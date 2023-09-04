%% Controle
clear; clc; close all;

%--------------------------------------------------------------------------
% Load some variables
%--------------------------------------------------------------------------
VARS   = load('./Results/Z-BinauralData_B.mat','sp_paramt','no_paramt');

%--------------------------------------------------------------------------
% Directories and files
%--------------------------------------------------------------------------
DIRIN  = './Results/';                            % Z-BinauralData_B file directory
DIROUT = ['./Results/',strrep(VARS.no_paramt.NAMENOISE,'.wav',''),'/AFOE/',num2str( VARS.no_paramt.AZN ),'/'];  % output file directory
NAMOUT = 'Result';

%--------------------------------------------------------------------------
% Adaptive Filter with Online Estimation Parameters
%--------------------------------------------------------------------------
beta  = 2e-2;   % MWF coefficient
gamma = 0; %era 3e-2
lambda_in      = 0.999;  % smothing coeffitient for input and speech, olny used when ( gamma ~= 0)
lambda_no      = 0.995;  % smothing coeffitient for noise ( gamma ~= 0)

%--------------------------------------------------------------------------
% Control Variables
%--------------------------------------------------------------------------
log               = 1;     % control the information on the screen (0 = on or 1 = off).
plot              = 0;     % plot (0  = on or 1 = off).
prior_est_time    = 6;     % time used to estimate the first estimatives of the correlation matrices (6 sec., default for my tests).
measurements      = 2;     % calculate measures: 0 = unable, 1 = only measure, 2 = measure and write audios
measurements_time = 5;     % time final segment of input signal used to measure de all the metrics in sec (4 sec, default for my tests).

af_paramt = struct( 'Beta', beta, 'Gamma', gamma,... % Coloquei -1 pra gamma so pra inicializar com alguma coisa
    'Lambda_in', lambda_in,  'Lambda_no', lambda_no,...
    'Log', log, 'Plot', plot, 'MeasurementsTime', measurements_time,...
    'PriorEstTime', prior_est_time, 'Measurements', measurements );

%--------------------------------------------------------------------------
% Input and output files
%--------------------------------------------------------------------------
NAMIN   = 'Z-BinauralData_B';   % input data filename

% files and load
FnamIn  = [ DIRIN, NAMIN ];     % input data
load( FnamIn );                 % load input data


%% Processamento

if( af_paramt.Log )
    disp('E: ONLINE FILTER');
end

% fft parameters ----------------------------------------------------------
M  =   256;         % length of the fft (48kHz and M=256 W=32 S=50)
W  =   128;         % length of the analysis window (128) [1]
S  =    50;         % overlap (percentual) (50%) [1]

% output variables --------------------------------------------------------

out_frq  = zeros(M,8);   % 1 - in_sp_left    5 - in_sp_right
pow_frq  = zeros(M,8);   % 2 - in_no_left    6 - in_no_right
itd_sp   = zeros(M,1);   % 3 - ou_sp_left    7 - ou_sp_right
itd_no   = zeros(M,1);   % 4 - ou_no_left    8 - ou_no_righ
ild_sp   = zeros(M,1);
ild_no   = zeros(M,1);

itd_sp_o = zeros(M,1);
itd_sp_i = zeros(M,1);
itd_no_o = zeros(M,1);
itd_no_i = zeros(M,1);

total_frames      = floor( tamsign/(W*(1-S/100)) ); 
pow_frq_inst      = zeros(M,8);
Delta_ild_no_inst = zeros( total_frames, 1 );
measurements = struct( 'out_frq', out_frq, 'pow_frq', pow_frq, 'itd_sp', itd_sp,...
                       'itd_no', itd_no, 'ild_sp', ild_sp, 'ild_no', ild_no,...
                       'itd_sp_o', itd_sp_o, 'itd_sp_i', itd_sp_i,...
                       'itd_no_o', itd_no_o,'itd_no_i', itd_no_i,...
                       'Time_dif_arrival_sp', 0, 'Time_dif_arrival_no', 0,...
                       'Snr_in_frt_left', 0, 'Snr_in_frt_right', 0,...
                       'Snr_ou_left', 0, 'Snr_ou_right', 0,...
                       'snr_in_left', 0, 'snr_in_right', 0,...
                       'snr_ou_left', 0, 'snr_ou_right', 0,...
                       'Delta_itd_sp', 0, 'Delta_itd_no', 0,...
                       'Delta_ild_sp', 0, 'Delta_ild_no', 0,...
                       'Pesq_in_left', 0, 'Pesq_in_right', 0,...
                       'Pesq_ou_left', 0, 'Pesq_ou_right', 0,...
                       'pow_frq_inst', pow_frq_inst,...
                       'ild_in_sp',0,'ild_ou_sp',0,...
                       'ild_in_no',0,'ild_ou_no',0);

%**************************************************************************
% Adaptive parameters
%**************************************************************************
af_coef        = zeros(2*numchan,M);            % adaptive filter coefficients
af_coef(1,:)   = 1;
af_coef(3/2*numchan+1,:) = 1;

%**************************************************************************
% Some vectors and matrices
%**************************************************************************

qL = [ 1 ; zeros(numchan-1,1) ];                        % left reference channel
qR = [ zeros(numchan/2,1) ; 1 ; zeros(numchan/2-1,1) ]; % right reference channel

%*************************************************************************
% Correlation matrices
%*************************************************************************

% weighted overlap-add memory allocation ------------------------------
buff_in_tmp_in = zeros(M,numchan);    % input buffer in the time domain
buff_in_frq_in = zeros(M,numchan);    % input in the frequency domain
a = zeros(M,numchan);    % input in the frequency domain
b = zeros(M,numchan);    % input in the frequency domain
buff_in_tmp_sp = zeros(M,numchan);    % input buffer in the time domain
buff_in_frq_sp = zeros(M,numchan);    % input in the frequency domain
buff_in_tmp_no = zeros(M,numchan);    % input buffer in the time domain
buff_in_frq_no = zeros(M,numchan);    % input in the frequency domain
buff_in_frq_no_ref  = zeros(M,2);     % noise signal in the reference microphone (front)
buff_in_frq_sp_ref  = zeros(M,2);     % speech signal in the reference microphone (front)
buff_out_tmp_in  = zeros(M,2);        % output buffer in the time domain (input signal)
buff_out_tmp_in_aux  = zeros(M,2);
c = zeros(M,2);
d = zeros(M,2);
buff_out_frq_in  = zeros(M,2);        % output buffer in the frequency domain  (input signal)
buff_out_tmp_sp  = zeros(M,2);        % output buffer in the time domain (speech signal)
buff_out_frq_sp  = zeros(M,2);        % output buffer in the frequency domain (speech signal)
buff_out_tmp_no  = zeros(M,2);        % output buffer in the time domain (noise signal)
buff_out_frq_no  = zeros(M,2);        % output buffer in the frequency domain (noise signal)
output_signal    = zeros(tamsign,2);  % output signal
output_speech    = zeros(tamsign,2);  % output speech
output_noise     = zeros(tamsign,2);  % output noise
buff_vad     = zeros(M,1);            % vad buffer
wm           = exp(1i*2*pi/M);        % complex exponential
bloco        = round(W*(1-S/100));    % innovation block
janela       = [zeros((M-W)/2,1);flip(sqrt(hann(W)));zeros((M-W)/2,1)];
count_frames = 0;
contv        = 0;
conty        = 0;
frameEndPriorEst  = round(af_paramt.PriorEstTime*FAM/(W*(1-S/100)))+1;
segment_length = af_paramt.MeasurementsTime*FAM;
segment        = tamsign - segment_length + 1;

% correlation matrices memory allocation ----------------------------------
Rx        = zeros(numchan,numchan,M);    % speech correlation matrix
Ry        = zeros(numchan,numchan,M);    % noisy signal correlation matrix
Rv        = zeros(numchan,numchan,M);    % noise correlation matrix

onlyOne = 0;
for sample = 1:tamsign
    
    buff_vad = [ speech_vad(sample) ; buff_vad(1:M-1) ];
    vad = sum(buff_vad); % vad's choice
    
    buff_in_tmp_in = [ buff_in_tmp_in(2:M,:) ; input_signal(sample,:) ];
    buff_in_tmp_sp = [ buff_in_tmp_sp(2:M,:) ; speech(sample,:) ];
    buff_in_tmp_no = [ buff_in_tmp_no(2:M,:) ; noise(sample,:) ];
    
    buff_out_tmp_in(:,1) = [ buff_out_tmp_in(2:M,1) ; 0 ];
    buff_out_tmp_in(:,2) = [ buff_out_tmp_in(2:M,2) ; 0 ];
    buff_out_tmp_sp(:,1) = [ buff_out_tmp_sp(2:M,1) ; 0 ];
    buff_out_tmp_sp(:,2) = [ buff_out_tmp_sp(2:M,2) ; 0 ];
    buff_out_tmp_no(:,1) = [ buff_out_tmp_no(2:M,1) ; 0 ];
    buff_out_tmp_no(:,2) = [ buff_out_tmp_no(2:M,2) ; 0 ];
    
    if( ~mod(sample,bloco) )
        
        % phase modification ----------------------------------------------
        count_frames = count_frames + 1;
        fact  = (((-1).^(0:M-1)).*(wm.^(-count_frames*bloco*(0:M-1))))';
        %if (count_frames <= 1000)
        for i = 1:256
                        if (abs(real(fact(i))) < 0.999999999 && abs(real(fact(i))) > 0.0000000001) 
                            disp(count_frames);
                            disp(fact(i));
                            disp('##############################');
                        else
                            if (abs(imag(fact(i))) < 0.999999999 && abs(imag(fact(i))) > 0.0000000001) 
                                disp(count_frames);
                                disp(fact(i));
                                disp('##############################');
                            end
                        end
              
        end
        % analysis --------------------------------------------------------
        for channel = 1 : numchan
            a(:, channel) = janela.*buff_in_tmp_in(:,channel);
            b(:, channel) = fft(janela.*buff_in_tmp_in(:,channel));
            buff_in_frq_in(:,channel) = fact.*b(:, channel);
           % buff_in_frq_in(:,channel) = fact.*fft(janela.*buff_in_tmp_in(:,channel));
            buff_in_frq_sp(:,channel) = fact.*fft(janela.*buff_in_tmp_sp(:,channel));
            buff_in_frq_no(:,channel) = fact.*fft(janela.*buff_in_tmp_no(:,channel));      
        end
        
        if( count_frames <= frameEndPriorEst ) 
            if( ~mod(sample,1000) && af_paramt.Log )
                disp(['   Prior estimation of corr. matrices: ', num2str(floor(100*sample/tamsign)),'%']);
            end
            
            % Correlation matrices prior estimative
            if( sum(sum( buff_in_frq_in )) ~= 0 )
            
                if( vad == 0 ) % noise frames
                    for bin = 1:M
                        Rv(:,:,bin) = Rv(:,:,bin) + (buff_in_frq_in(bin,:).') * conj(buff_in_frq_in(bin,:));
                    end
                    contv = contv + 1; % number of noise frames
                    
                else % speech + noise frames
                    for bin = 1:M
                        Ry(:,:,bin) = Ry(:,:,bin) + (buff_in_frq_in(bin,:).') * conj(buff_in_frq_in(bin,:));
                    end
                    conty = conty + 1; % number of speech + noise frames
                end
                
            end
            
            % Filtering: trivial filter
            for bin = 1:M
                buff_out_frq_in(bin,1) = qL'*(buff_in_frq_in(bin,:).');
                buff_out_frq_in(bin,2) = qR'*(buff_in_frq_in(bin,:).');
                
                buff_out_frq_sp(bin,1) = qL'*(buff_in_frq_sp(bin,:).');
                buff_out_frq_sp(bin,2) = qR'*(buff_in_frq_sp(bin,:).');
                
                buff_out_frq_no(bin,1) = qL'*(buff_in_frq_no(bin,:).');
                buff_out_frq_no(bin,2) = qR'*(buff_in_frq_no(bin,:).');
            end
            
            
        else
            
            % Normalization: performed only once
            if( onlyOne < 1 )
                Ry = Ry / (conty-1);
                Rv = Rv / (contv-1);

                Rx = Ry - Rv;
                
                onlyOne = onlyOne + 1; % make if clause false
            end
            
            %--------------------------------------------------------------
            % Online filtering and Quality measurements
            %--------------------------------------------------------------
            if( ~mod(sample,1000) && af_paramt.Log )
                disp([ '   Filtering: ', num2str( floor( 100*sample/tamsign )), '%' ]);
            end
            
            
            if( ( af_paramt.Lambda_in >= 0 || af_paramt.af_paramt.Lambda_no >= 0 ) )

                % Update correlation matrices
                if( vad == 0  ) % noise
                    for bin = 1:M
                        Rv(:,:,bin) = af_paramt.Lambda_no * Rv(:,:,bin) + ( 1 - af_paramt.Lambda_no ) *( buff_in_frq_in(bin,:).' * conj(buff_in_frq_in(bin,:)));
                    end
                else % speech + noise
                    
                    for bin = 1:M 
                        Ry(:,:,bin) = af_paramt.Lambda_in * Ry(:,:,bin) + ( 1 - af_paramt.Lambda_in ) *( buff_in_frq_in(bin,:).' * conj(buff_in_frq_in(bin,:)) );
                        Rx(:,:,bin) = Ry(:,:,bin) - Rv(:,:,bin);
                    end
                end
                
                % Update adaptive coefficients
                af_coef = Binaural_D_Search_MWFILD( af_paramt, af_coef, Rx, Ry, Rv, qL, qR, M );
            end
            
            % Filtering: output of the multichannel Wiener filter
            for bin = 1:M
                
                % processed input signal frame
                buff_out_frq_in(bin,1) = af_coef(1:numchan,bin)' * (buff_in_frq_in(bin,:).');
                buff_out_frq_in(bin,2) = af_coef(numchan+1:2*numchan,bin)' * (buff_in_frq_in(bin,:).');
                
                % processed speech signal frame
                buff_out_frq_sp(bin,1) = af_coef(1:numchan,bin)' * (buff_in_frq_sp(bin,:).');
                buff_out_frq_sp(bin,2) = af_coef(numchan+1:2*numchan,bin)' * (buff_in_frq_sp(bin,:).');
                
                % processed noise signal frame
                buff_out_frq_no(bin,1) = af_coef(1:numchan,bin)'* (buff_in_frq_no(bin,:).');
                buff_out_frq_no(bin,2) = af_coef(numchan+1:2*numchan,bin)'* (buff_in_frq_no(bin,:).');
                
                % speech in reference microphone ( front ) 
                buff_in_frq_sp_ref(bin,1) = qL'*(buff_in_frq_sp(bin,:).');
                buff_in_frq_sp_ref(bin,2) = qR'*(buff_in_frq_sp(bin,:).');
                
                % noise in reference microphone ( front )
                buff_in_frq_no_ref(bin,1) = qL'*(buff_in_frq_no(bin,:).');
                buff_in_frq_no_ref(bin,2) = qR'*(buff_in_frq_no(bin,:).');
            end

            % Objective measurements
            if( ( sample >= segment ) && ( sample <= fim ) )
                measurements = Binaural_F_QualityMeasurements_AFOE( buff_in_frq_sp_ref, buff_in_frq_no_ref, buff_out_frq_sp, buff_out_frq_no, measurements, M, FAM, count_frames );
            end
        end
        
        % Synthesis
        c(:, 1) = conj(fact).*buff_out_frq_in(:,1);
        d(:, 1) = ifft(conj(fact).*buff_out_frq_in(:,1));
        buff_out_tmp_in_aux(:, 1) = janela .* ifft(conj(fact).*buff_out_frq_in(:,1));
        buff_out_tmp_in(:,1) = buff_out_tmp_in(:,1) + buff_out_tmp_in_aux(:, 1);
       % buff_out_tmp_in(:,1) = buff_out_tmp_in(:,1) + janela .* ifft(conj(fact).*buff_out_frq_in(:,1));
        c(:, 2) = conj(fact).*buff_out_frq_in(:,2);
        d(:, 2) = ifft(conj(fact).*buff_out_frq_in(:,2));
        buff_out_tmp_in_aux(:, 2) = janela .* ifft(conj(fact).*buff_out_frq_in(:,2));
        buff_out_tmp_in(:,2) = buff_out_tmp_in(:,2) + buff_out_tmp_in_aux(:, 2);
        %buff_out_tmp_in(:,2) = buff_out_tmp_in(:,2) + janela .* ifft(conj(fact).*buff_out_frq_in(:,2));

        buff_out_tmp_sp(:,1) = buff_out_tmp_sp(:,1) + janela .* ifft(conj(fact).*buff_out_frq_sp(:,1));
        buff_out_tmp_sp(:,2) = buff_out_tmp_sp(:,2) + janela .* ifft(conj(fact).*buff_out_frq_sp(:,2));

        buff_out_tmp_no(:,1) = buff_out_tmp_no(:,1) + janela .* ifft(conj(fact).*buff_out_frq_no(:,1));
        buff_out_tmp_no(:,2) = buff_out_tmp_no(:,2) + janela .* ifft(conj(fact).*buff_out_frq_no(:,2));

        % Overlap-add
        output_signal(sample-bloco+1:sample,1) = real(buff_out_tmp_in(1:bloco,1));
        output_signal(sample-bloco+1:sample,2) = real(buff_out_tmp_in(1:bloco,2));

        output_speech(sample-bloco+1:sample,1) = real(buff_out_tmp_sp(1:bloco,1));
        output_speech(sample-bloco+1:sample,2) = real(buff_out_tmp_sp(1:bloco,2));

        output_noise(sample-bloco+1:sample,1) = real(buff_out_tmp_no(1:bloco,1));
        output_noise(sample-bloco+1:sample,2) = real(buff_out_tmp_no(1:bloco,2));
    end
    if (count_frames == 2)
        return;
    end
end

disp('   Measurement');
                
% output files
NAMUNIN = 'B1-unprocessed_input';    % unprocessed input filename 
NAMPRIN = 'B2-processed_input';      % processed input filename 
NAMUNSP = 'B3-unprocessed_speech';   % unprocessed speech filename
NAMPRSP = 'B4-processed_speech';     % processed speech filename
NAMUNNO = 'B5-unprocessed_noise';    % unprocesssed noise filename
NAMPRNO = 'B6-processed_noise';      % processsed noise filename
NAMCPNO = 'B7-comparison_noise';     % comparison noise filename

% output files: full signals
Sx = [ 'S', num2str(sp_paramt.AZS) ];
Nx = [ 'N', num2str(no_paramt.AZN) ];
GAMMASTR = [ num2str(af_paramt.Gamma,'%10.1e') ];

FnamUnIn = [DIROUT,NAMUNIN,'_',Sx,Nx,'.wav'];
FnamUnInSeg = [DIROUT,NAMUNIN,'_',Sx,Nx,'_Seg.wav'];
FnamUnSpSeg = [DIROUT,NAMUNSP,'_',Sx,'_Seg.wav'];
FnamUnNoSeg = [DIROUT,NAMUNNO,'_',Nx,'_Seg.wav'];

if( af_paramt.Gamma == 0 ) % name when GAMMA = 0
    FnamPrIn = [DIROUT,NAMPRIN,'_MWF_',Sx,Nx,'.wav'];   % processed input 
    FnamPrInSeg = [DIROUT,NAMPRIN,'_MWF_',Sx,Nx,'_Seg.wav'];   % processed input 
    FnamPrSpSeg = [DIROUT,NAMPRSP,'_MWF_',Sx,'_Seg.wav'];      % processed speech
    FnamPrNoSeg = [DIROUT,NAMPRNO,'_MWF_',Nx,'_Seg.wav'];      % processsed noise
    FnamCpNoSeg = [DIROUT,NAMCPNO,'_MWF_',Nx,'_Seg.wav'];      % comparison noise
else
    FnamPrIn = [DIROUT,NAMPRIN,'_',Sx,Nx,'_MWFILD_',GAMMASTR,'.wav'];   % processed input 
    FnamPrInSeg = [DIROUT,NAMPRIN,'_',Sx,Nx,'_MWFILD_',GAMMASTR,'_Seg.wav'];   % processed input 
    FnamPrSpSeg = [DIROUT,NAMPRSP,'_',Sx,'_MWFILD_',GAMMASTR,'_Seg.wav'];      % processed speech
    FnamPrNoSeg = [DIROUT,NAMPRNO,'_',Nx,'_MWFILD_',GAMMASTR,'_Seg.wav'];      % processsed noise
    FnamCpNoSeg = [DIROUT,NAMCPNO,'_',Nx,'_MWFILD_',GAMMASTR,'_Seg.wav'];      % comparison noise
end
output_signal_seg = output_signal( segment:tamsign, : );
output_speech_seg = output_speech( segment:tamsign, : );
output_noise_seg  = output_noise( segment:tamsign, : );
input_signal_seg  = input_signal( segment:tamsign, : );
speech_seg        = speech( segment:tamsign,: );
noise_seg         = noise( segment:tamsign, : );

% output noise with same power as input noise -----------------------------
comp_noise_seg = zeros(length(noise_seg),2);
aux_left  = std( noise_seg( :, 1 ));
aux_right = std( noise_seg( :, numchan/2 + 1 ));
if( aux_left > aux_right )
    aux =  aux_left / std( output_noise_seg( :, 1 ));
else
    aux =  aux_right / std( output_noise_seg( :, 2 ));
end
comp_noise_seg(:,1) = output_noise_seg(:,1) * aux;
comp_noise_seg(:,2) = output_noise_seg(:,2) * aux;

audiowrite( FnamPrInSeg, output_signal_seg, FAM  );
audiowrite( FnamPrSpSeg, output_speech_seg, FAM );
audiowrite( FnamPrNoSeg, output_noise_seg, FAM );
audiowrite( FnamCpNoSeg, comp_noise_seg, FAM );
audiowrite( FnamUnInSeg, [input_signal_seg(:,1) input_signal_seg(:,numchan/2+1)], FAM );
audiowrite( FnamUnSpSeg, [speech_seg(:,1) speech_seg(:,numchan/2+1)], FAM );
audiowrite( FnamUnNoSeg, [noise_seg(:,1) noise_seg(:,numchan/2+1)], FAM );

audiowrite( FnamPrIn, output_signal, FAM  );
audiowrite( FnamUnIn, [input_signal(:,1) input_signal(:,numchan/2+1)], FAM );


%**************************************************************************
% Quality Measurement
%**************************************************************************

cd ./Progs
measurements.Pesq_in_left  = pesqbin(speech_seg(:,1),input_signal_seg(:,1),FAM,'wb');
measurements.Pesq_in_right = pesqbin(speech_seg(:,4),input_signal_seg(:,4),FAM,'wb');
measurements.Pesq_ou_left  = pesqbin(speech_seg(:,1),output_signal_seg(:,1),FAM,'wb');
measurements.Pesq_ou_right = pesqbin(speech_seg(:,4),output_signal_seg(:,2),FAM,'wb');
cd ..
measurements.Delta_Pesq_left  = measurements.Pesq_ou_left  - measurements.Pesq_in_left;
measurements.Delta_Pesq_right = measurements.Pesq_ou_right - measurements.Pesq_in_right;

%**************************************************************************
% Time difference of arrival [3]
%**************************************************************************
cd ./Progs
time_sp_in = sigalign(speech_seg(:,1),speech_seg(:,4));
time_no_in = sigalign(noise_seg(:,1),noise_seg(:,4));
time_sp_ou = sigalign(output_speech_seg(:,1),output_speech_seg(:,2));
time_no_ou = sigalign(output_noise_seg(:,1),output_noise_seg(:,2));
cd ..
measurements.Time_dif_arrival_sp = time_sp_ou - time_sp_in;
measurements.Time_dif_arrival_no = time_no_ou - time_no_in;

%**************************************************************************
% Input and output SNR from left and right channels
%**************************************************************************

% original and processed SNR ----------------------------------------------
measurements.Snr_in_frt_left  = 10*log10(sum(measurements.pow_frq(:,1))/sum(measurements.pow_frq(:,2)));
measurements.Snr_in_frt_right = 10*log10(sum(measurements.pow_frq(:,5))/sum(measurements.pow_frq(:,6)));
measurements.Snr_ou_left      = 10*log10(sum(measurements.pow_frq(:,3))/sum(measurements.pow_frq(:,4)));
measurements.Snr_ou_right     = 10*log10(sum(measurements.pow_frq(:,7))/sum(measurements.pow_frq(:,8)));

%**************************************************************************
% Interaural time difference [9]
%**************************************************************************

maxfrq = floor(1500*M/FAM+1);
measurements.Delta_itd_sp = (1/pi) * sum( abs( unwrap(angle(measurements.itd_sp_o(1:maxfrq))) - unwrap(angle(measurements.itd_sp_i(1:maxfrq))) ) ); %unnwrap(angle()) era phase()
measurements.Delta_itd_no = (1/pi) * sum( abs( unwrap(angle(measurements.itd_no_o(1:maxfrq))) - unwrap(angle(measurements.itd_no_i(1:maxfrq))) ) ); %unnwrap(angle()) era phase()

%**************************************************************************
% Interaural level difference [9]
%**************************************************************************

minfrq = floor(0*M/FAM+1);
aux_sp_in = (measurements.pow_frq(:,1)./measurements.pow_frq(:,5));
aux_no_in = (measurements.pow_frq(:,2)./measurements.pow_frq(:,6));

aux_sp_ou = (measurements.pow_frq(:,3)./measurements.pow_frq(:,7));
aux_no_ou = (measurements.pow_frq(:,4)./measurements.pow_frq(:,8));

ILD_in_sp = 10*log10(aux_sp_in(minfrq:M/2));
ILD_ou_sp = 10*log10(aux_sp_ou(minfrq:M/2));

ILD_in_no = 10*log10(aux_no_in(minfrq:M/2));
ILD_ou_no = 10*log10(aux_no_ou(minfrq:M/2));

measurements.ild_in_sp = (1/(M/2+1-minfrq)) * sum(ILD_in_sp);
measurements.ild_ou_sp = (1/(M/2+1-minfrq)) * sum(ILD_ou_sp);

measurements.ild_in_no = (1/(M/2+1-minfrq)) * sum(ILD_in_no);
measurements.ild_ou_no = (1/(M/2+1-minfrq)) * sum(ILD_ou_no);

measurements.Delta_ild_sp = abs(measurements.ild_ou_sp - measurements.ild_in_sp);
measurements.Delta_ild_no = abs(measurements.ild_ou_no - measurements.ild_in_no);


%% Metricas
measures = zeros(14,1);
measures( 1, 1 )  = measurements.Pesq_ou_left;
measures( 2, 1 )  = measurements.Pesq_ou_right;
measures( 3, 1 )  = measurements.Delta_Pesq_left;
measures( 4, 1 )  = measurements.Delta_Pesq_right;
measures( 5, 1 )  = measurements.Snr_ou_left;
measures( 6, 1 )  = measurements.Snr_ou_right;
measures( 7, 1 ) = measurements.Delta_itd_sp;
measures( 8, 1 ) = measurements.Delta_itd_no;
measures( 9, 1 ) = measurements.ild_in_sp;
measures( 10, 1 ) = measurements.ild_ou_sp;
measures( 11, 1 ) = measurements.ild_in_no;
measures( 12, 1 ) = measurements.ild_ou_no;
measures( 13, 1 ) = measurements.Delta_ild_sp;
measures( 14, 1 ) = measurements.Delta_ild_no;


FnamOutSeg  = [DIROUT, 'Result_AFOE_', 'S', num2str(sp_paramt.AZS),'N', num2str(no_paramt.AZN),'_Seg'];

save( FnamOutSeg, 'measures' );

disp( measures )
beep


%% Other functions used above


function measurements = Binaural_F_QualityMeasurements_AFOE( frame_in_sp_ref, frame_in_no_ref, frame_out_sp, frame_out_no, measurements, M, FAM, count_frames )

% speech in front microphones -----------------------------------------
measurements.out_frq(:,1) = frame_in_sp_ref(:,1);
measurements.out_frq(:,5) = frame_in_sp_ref(:,2);

% noise in front microphones ------------------------------------------
measurements.out_frq(:,2) = frame_in_no_ref(:,1);
measurements.out_frq(:,6) = frame_in_no_ref(:,2);

% processed speech ----------------------------------------------------
measurements.out_frq(:,3) = frame_out_sp(:,1);
measurements.out_frq(:,7) = frame_out_sp(:,2);

% processed noise -----------------------------------------------------
measurements.out_frq(:,4) = frame_out_no(:,1);
measurements.out_frq(:,8) = frame_out_no(:,2);

% input (front) and output ------------------------------------------------
for bin = 1 : M
    
    % itd -----------------------------------------------------------------
    measurements.itd_sp(bin) = measurements.itd_sp(bin) + ...
        ( unwrap(angle(measurements.out_frq(bin,3) * ...                %unwrap(angle()) era phase()
        conj(measurements.out_frq(bin,7)) * measurements.out_frq(bin,5) * ...    
        conj(measurements.out_frq(bin,1)) )) )^2;

    measurements.itd_no(bin)  = measurements.itd_no(bin) + ...
        ( unwrap(angle(measurements.out_frq(bin,4) * ...                %unwrap(angle()) era phase()
        conj(measurements.out_frq(bin,8)) * measurements.out_frq(bin,6) * ...
        conj(measurements.out_frq(bin,2)) )) )^2;

    measurements.itd_sp_o(bin) = measurements.itd_sp_o(bin) + ...
        measurements.out_frq(bin,3) * conj(measurements.out_frq(bin,7));

    measurements.itd_sp_i(bin) = measurements.itd_sp_i(bin) + ...
        measurements.out_frq(bin,1) * conj(measurements.out_frq(bin,5));

    measurements.itd_no_o(bin) = measurements.itd_no_o(bin) + ...
        measurements.out_frq(bin,4) * conj(measurements.out_frq(bin,8));

    measurements.itd_no_i(bin) = measurements.itd_no_i(bin) + ...
        measurements.out_frq(bin,2) * conj(measurements.out_frq(bin,6));

    % ild -----------------------------------------------------
    measurements.ild_sp(bin) = measurements.ild_sp(bin) + ...
        ( log10(abs(measurements.out_frq(bin,3) * ...
        conj(measurements.out_frq(bin,5))) / abs(measurements.out_frq(bin,7) * ...
        conj(measurements.out_frq(bin,1)))) )^2;

    measurements.ild_no(bin) = measurements.ild_no(bin) + ...
        ( log10(abs(measurements.out_frq(bin,4) * ...
        conj(measurements.out_frq(bin,6))) / abs(measurements.out_frq(bin,8) * ...
        conj(measurements.out_frq(bin,2)))) )^2;
end

% power -------------------------------------------------------------------
measurements.pow_frq = measurements.pow_frq + abs(measurements.out_frq).^2;

% Instantaneous
measurements.pow_frq_inst = abs( measurements.out_frq ).^2; % average power           

% Interaural Level Difference
minfrq = floor(0*M/FAM+1);
aux_no = ( measurements.pow_frq_inst(:,4).*measurements.pow_frq_inst(:,6))./(measurements.pow_frq_inst(:,2).*measurements.pow_frq_inst(:,8));
ild_no = abs(10*log10(aux_no(minfrq:M/2)));
measurements.Delta_ild_no_inst( count_frames ) = (1/( M/2+1-minfrq )) * sum(ild_no);
end

