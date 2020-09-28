%% Super Lorentzian lineshape including second moment, and also interpolated.
% Return LUT interpolated from -30kHz to +30KHz.
% Units of Gs. Shaihan Malik and Daniel West 2019.

function [G,D,ff] = SuperLorentzian_LSint(T2b,fsample)

% Define frequency range.
n=512;
ff = linspace(-30e3,30e3,n);

% Compute G for this range.
G = zeros([n 1]);
for ii=1:n
    G(ii) = SL(ff(ii));
end

% Interpolate.
po = find(abs(ff)<1.5e3); % Points to interpolate.
pu = find((abs(ff)>1.5e3)&(abs(ff)<2e3)); % Points to use.

Gi = spline(ff(pu),G(pu),ff(po));
G(po) = Gi;

D = 1/(sqrt(15)*T2b);

% If  user requested samples, sample the correct values.
if nargin==2
    G = interp1(ff,G,fsample);
end

    function gg = SL(f)
        th = linspace(0,pi/2,500);
        dth = th(2)-th(1);
        g = sin(th).*sqrt(2/pi).*T2b./(abs(3*cos(th).^2-1));
        g = g .* exp(-2*(2*pi*f*T2b./abs(3*cos(th).^2-1)).^2);
        gg = dth*sum(g);
    end

end


