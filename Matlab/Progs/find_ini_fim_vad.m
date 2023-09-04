clear; clc;
load('C:\Users\User11\Dropbox\Ufsc_dropbox\Simulacoes\ProgEusipco2014\Results\Z-BinauralData_B_n60_p60_n60_p60_3s_1s.mat')

% search for the initial nonzero vad sample -------------------------------
flag    = 1;    
cont    = 0;
for amostra = 1:tamsign
   
    if( (speech_vad( amostra ) > 0) && ( flag == 1 ) )
        cont        = cont + 1;
        ini( cont ) =  amostra ;
        flag        = 0;
    elseif( (speech_vad( amostra ) == 0) && ( flag == 0 ) )
        flag = 1;
    end
end

% search for the final nozero vad sample ----------------------------------
flag    = 1;
cont    = 0;
for amostra = 1:tamsign
    if( (speech_vad( amostra ) > 0) && ( flag == 1 ) )
        flag = 0;
    elseif( (speech_vad( amostra ) == 0) && ( flag == 0 ) )
        cont        = cont + 1;
        fim( cont ) = amostra;
        flag        = 1;
    end
end

ini(end-3)
fim(end)

