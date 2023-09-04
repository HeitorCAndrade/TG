%**************************************************************************
% Binaural_D_Search_AF_v2_test
%
%
%
%
% Author   : Diego Marques do Carmo
% Created  : 13/02/2015
%**************************************************************************

function w = Binaural_D_Search_MWFILD( af_paramt, initial_coeffs, RX, RY, RV, qL, qR, M )

w  = initial_coeffs;
CH = length(qL);
q  = [qL; qR];              % vector of references microphone
Zm = zeros(CH,CH);
Imm = eye(2*CH,2*CH);
Im  = eye(CH,CH);
for bin = 1:M/2+1
    
    Ry = RY(:,:,bin);
    Rx = RX(:,:,bin);
    Rv = RV(:,:,bin);
    
    Rxx = [Rx, Zm; Zm, Rx];
    Ryy = [Ry, Zm; Zm, Ry];
    
    % Constant parts in the JMWF function
    v1 = (af_paramt.Beta*Rxx*q);
    M1 = Imm - (af_paramt.Beta*Ryy);
    
    % MWF part - this part has the previous coefficients
    w_MWF = M1*w(:,bin) + v1;
    if( af_paramt.Gamma == 0)
        
        w(:,bin) = w_MWF;
    else
        R1  = [ Rv, Zm; Zm, Zm];
        R2  = [ Zm, Zm; Zm, Rv];
        Rvv = R1 + R2;
    
        ccL = qL'*Rv*qL;
        ccR = qR'*Rv*qR;
    
        % ILD part
        pLn = w(:,bin)'*R1*w(:,bin);
        pRn = w(:,bin)'*R2*w(:,bin);
        dn  = pLn*ccR + pRn*ccL;
        a_1 = ccR*ccL^2*pRn^2;
        a_2 = ccL*ccR^2*pLn^2;
        b_1 = ccL*ccR^2*pLn*pRn;
        b_2 = ccR*ccL^2*pLn*pRn;
        Cn  = [(a_1-b_1)*Im Zm; Zm (a_2-b_2)*Im];
        
        w_ILD = (Cn*Rvv*w(:,bin))/dn^3;
        
        w(:,bin) = w_MWF + af_paramt.Gamma*w_ILD;
    end
    
    if (bin > 1) && (bin <= M/2)
        w(:,M+2-bin) = conj(w(:,bin));
    end
end