%**************************************************************************
% Binaural_D_Search_AF_v2_test
%
%
%
%
% Author   : Diego Marques do Carmo
% Created  : 13/02/2015
%**************************************************************************

function w = Binaural_D_Online_MWF( af_paramt, initial_coeffs, Rx, Ry, qL, qR )

w = initial_coeffs;

q = [qL; qR];              % vector of references microphone

Zm  = zeros(length(qL),length(qR));

Rxx = [ Rx, Zm;
    Zm, Rx];

Ryy = [ Ry, Zm;
    Zm, Ry];

I   = eye(size(Ryy));

% Constant parts in the JMWF function
v1 = (af_paramt.Beta*Rxx*q);
M1 = I - (af_paramt.Beta*Ryy);
w = M1*w + v1;

end