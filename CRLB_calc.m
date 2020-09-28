function [CRLB,FIM] = CRLB_calc(sig_diff)
%   function CRLB = CRLB_calc(sig_diff)
% Inputs:
%   sig_diff is a matrix [Npar x Nmeas]
% Outputs:
%   CRLB is a column vector with variances

[Npar, Nmeas] = size(sig_diff);

FIM = zeros(Npar);

for nn=1:Nmeas
    
    for jj=1:Npar
        for kk=1:Npar
            FIM(jj,kk) = FIM(jj,kk) + sig_diff(jj,nn)*sig_diff(kk,nn);
        end
    end
    
end

CRLB = diag(pinv(FIM));  

end