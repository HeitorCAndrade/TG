%--------------------------------------------------------------------------
% Simulation_AFOE.m
%
% M.H. Costa, P.A. Naylor (2013) Interaural level difference
% preservation in the multichannel Wiener fielter for binaural hearing
% aid applications, Eusipco 2014.
%
% Author    : Diego Marques
%--------------------------------------------------------------------------

clear
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
FnamOutSeg  = [DIROUT, NAMOUT, '_AFOE_', 'S', num2str(VARS.sp_paramt.AZS),'N', num2str(VARS.no_paramt.AZN),'_Seg'];

%--------------------------------------------------------------------------
% Adaptive Filter with Online Estimation Parameters
%--------------------------------------------------------------------------
beta  = 2e-2;   % MWF coefficient
lambda_in      = 0.999;  % smothing coeffitient for input and speech, olny used when ( gamma ~= 0)
lambda_no      = 0.994;  % smothing coeffitient for noise ( gamma ~= 0)

%--------------------------------------------------------------------------
% Control Variables
%--------------------------------------------------------------------------
log               = 1;     % control the information on the screen (0 = on or 1 = off).
prior_est_time    = 6;     % time used to estimate the first estimatives of the correlation matrices (6 sec., default for my tests).
measurements_time = 4;     % time final segment of input signal used to measure de all the metrics in sec (4 sec, default for my tests).

af_paramt = struct( 'Beta', beta,...
    'Lambda_in', lambda_in, 'Lambda_no', lambda_no,...
    'Log', log, 'MeasurementsTime', measurements_time,...
    'PriorEstTime', prior_est_time );

%----------------------------------
% Processing
%----------------------------------

if( af_paramt.Log )
    disp('E: ONLINE FILTER');
end

%--------------------------------------------------------------------------
% Input and output files
%--------------------------------------------------------------------------
NAMIN   = 'Z-BinauralData_B';   % input data filename
NAMOUT  = ['Z-BinauralData_F_MWF_Seg.mat' ];  % output data filename

% files and load
FnamIn  = [ DIRIN, NAMIN ];     % input data
FnamOut = [ DIROUT, NAMOUT ];   % output data
load( FnamIn );                 % load input data

% fft parameters ----------------------------------------------------------
M  =   256;         % length of the fft
W  =   128;         % length of the analysis window (128) [1]
S  =    50;         % overlap (percentual) (50%) [1]

                   
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
buff_in_tmp_sp = zeros(M,numchan);    % input buffer in the time domain
buff_in_frq_sp = zeros(M,numchan);    % input in the frequency domain
buff_in_tmp_no = zeros(M,numchan);    % input buffer in the time domain
buff_in_frq_no = zeros(M,numchan);    % input in the frequency domain
buff_in_frq_no_ref  = zeros(M,2);     % noise signal in the reference microphone (front)
buff_in_frq_sp_ref  = zeros(M,2);     % speech signal in the reference microphone (front)
buff_out_tmp_in  = zeros(M,2);        % output buffer in the time domain (input signal)
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
janela       = [zeros((M-W)/2,1);flipdim(sqrt(hann(W)),2);zeros((M-W)/2,1)];
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
Rv_true   = zeros(numchan,numchan,M);    % ideal noise correlation matrix

onlyOne = 0;
for sample = 1:tamsign
    
    buff_vad = [ speech_vad(sample) ; buff_vad(1:M-1) ];
    vad = sum(buff_vad); % vad's choice
    
    buff_in_tmp_in = [ buff_in_tmp_in(2:M,:) ; input_signal(sample,:) ];
    
    buff_out_tmp_in(:,1) = [ buff_out_tmp_in(2:M,1) ; 0 ];
    buff_out_tmp_in(:,2) = [ buff_out_tmp_in(2:M,2) ; 0 ];
    
    if( ~mod(sample,bloco) )
        
        % phase modification ----------------------------------------------
        count_frames = count_frames + 1;
        fact  = (((-1).^(0:M-1)).*(wm.^(-count_frames*bloco*(0:M-1))))';
        
        % analysis --------------------------------------------------------
        for channel = 1 : numchan
            buff_in_frq_in(:,channel) = fact.*fft(janela.*buff_in_tmp_in(:,channel));
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
                        Rv(:,:,bin)      = af_paramt.Lambda_no * Rv(:,:,bin)      + ( 1 - af_paramt.Lambda_no ) *( buff_in_frq_in(bin,:).' * conj(buff_in_frq_in(bin,:)));
                    end
                else % speech + noise
                    
                    for bin = 1:M 
                        Ry(:,:,bin) = af_paramt.Lambda_in * Ry(:,:,bin) + ( 1 - af_paramt.Lambda_in ) *( buff_in_frq_in(bin,:).' * conj(buff_in_frq_in(bin,:)) );
                        Rx(:,:,bin) = Ry(:,:,bin) - Rv(:,:,bin);
                    end
                end
                
                % Update adaptive coefficients                
                for bin = 1:M/2+1
                    af_coef(:,bin) = Binaural_Online_MWF( af_paramt, af_coef(:,bin), Rx(:,:,bin), Ry(:,:,bin), qL, qR );
                    
                    % negative frequencies
                    if (bin > 1) && (bin <= M/2)
                        af_coef(:,M+2-bin) = conj(af_coef(:,bin));
                    end
                end
            end
            
            % Filtering: output of the multichannel Wiener filter
            for bin = 1:M
                
                % processed input signal frame
                buff_out_frq_in(bin,1) = af_coef(1:numchan,bin)' * (buff_in_frq_in(bin,:).');
                buff_out_frq_in(bin,2) = af_coef(numchan+1:2*numchan,bin)' * (buff_in_frq_in(bin,:).');
            end
        end
        
        % Synthesis
        buff_out_tmp_in(:,1) = buff_out_tmp_in(:,1) + janela .* ifft(conj(fact).*buff_out_frq_in(:,1));
        buff_out_tmp_in(:,2) = buff_out_tmp_in(:,2) + janela .* ifft(conj(fact).*buff_out_frq_in(:,2));

        % Overlap-add
        output_signal(sample-bloco+1:sample,1) = real(buff_out_tmp_in(1:bloco,1));
        output_signal(sample-bloco+1:sample,2) = real(buff_out_tmp_in(1:bloco,2));
    end
end


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

FnamUnInSeg = [DIROUT,NAMUNIN,'_',Sx,Nx,'_Seg.wav'];
FnamPrInSeg = [DIROUT,NAMPRIN,'_MWF_',Sx,Nx,'_Seg.wav'];   % processed input

FnamUnIn = [DIROUT,NAMUNIN,'_',Sx,Nx,'.wav'];
FnamPrIn = [DIROUT,NAMPRIN,'_MWF_',Sx,Nx,'.wav'];   % processed input
                
% write audios: last segment from the whole signal
output_signal_seg = output_signal( segment:tamsign, : );
input_signal_seg  = input_signal( segment:tamsign, : );

audiowrite( FnamPrInSeg, output_signal_seg, FAM  );
audiowrite( FnamUnInSeg, [input_signal_seg(:,1) input_signal_seg(:,numchan/2+1)], FAM );

audiowrite( FnamPrIn, output_signal, FAM  );
audiowrite( FnamUnIn, [input_signal(:,1) input_signal(:,numchan/2+1)], FAM );
