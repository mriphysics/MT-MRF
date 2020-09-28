function [pulses,pulseBB,pulseMF,pulse_weights,b1_max_pulse] = gen_CSMT_pulse_Diffnp(flips,dur,TR,b1rms_total,delta,nband,M,S,varargin)
%%% Function generates a set of pulses for each of the flips, with the
%%% specified B1rms for all pulses, and the specified duration. Uses a
%%% Gaussian shape. Adapted from Shaihan Malik 12-2-2018.

gausswin_sigma = 3; % Shape of basic pulse.
dt = 6.4e-6;

if ~isempty(varargin)
    for ii=1:length(varargin)
        if strcmpi(varargin{ii},'sigma')
            gausswin_sigma = varargin{ii+1};
        end
        
        if strcmpi(varargin{ii},'dt')
            dt = varargin{ii+1};
        end
        
    end
end

nflip = length(flips);

%% Set-up time domain.

nt = ceil(dur/dt);
dur = dt*nt;

%% Basic shape.

pulse = gausswin(nt,gausswin_sigma);
pulse = pulse - min(pulse);
pulse = pulse / max(pulse);
trms = sum(pulse.^2)*dt/max(pulse)^2; % [s]
gam = 267.5221; % [rad/s/uT]

teff =  (gam*sum(pulse)*dt)/(gam*max(pulse))/dur; % teff of shape - normalised to duration.

%% Compute ratio of off to on resonance power. 
% Maximum flip angle sets the max B1 which is the constraint here.

b1_max_ONres = max(flips)/(gam*teff*dur); % [uT]
b1_rms_ONres = sqrt(b1_max_ONres^2 * trms / TR); % RMS over the TR.

% Compute k factor.
switch nband
    case 2
        k = ((M+S)/(4*M))*((b1rms_total/b1_rms_ONres)^2-1);
    case 3
        k = ((M+S)/(2*M))*((b1rms_total/b1_rms_ONres)^2-1);
end

% Modulation function maximum.
mfmax = (2*sqrt(k)+1);

% Compute max pulse amplitude.
b1_max_pulse = mfmax * b1_max_ONres;


%% Make pulses - repeat the above calculations again for each band.

pulses = {};
tt = dt*(1:nt);

for jj=1:nflip
    flip_current = flips(jj);
      
    b1_max_ONres = (flip_current)/(gam*teff*dur); % [uT]
    b1_rms_ONres = sqrt(b1_max_ONres^2 * trms / TR); % RMS over the TR.
    
    switch nband
        
        case 2

            % Compute k factor.
            k = ((M+S)/(4*M))*((b1rms_total/b1_rms_ONres)^2-1);
            
            % Modulation function maximum.
            mfmax = (2*sqrt(k)+1);
            
            % Compute max pulse amplitude.
            b1_max_pulse = mfmax * b1_max_ONres;

            % Make modulated pulse.
            mf = 2*sqrt(k)*(cos(2*pi*delta*tt)+1i*sin(2*pi*delta*tt)) + 1;
            pulses{jj} = pulse(:) .* mf(:) * b1_max_ONres;
            
        case 3
            
            % Compute k factor.
            k = ((M+S)/(2*M))*((b1rms_total/b1_rms_ONres)^2-1);
            
            % Modulation function maximum.
            mfmax = (2*sqrt(k)+1);
            
            % Compute max pulse amplitude.
            b1_max_pulse = mfmax * b1_max_ONres;

            % Make modulated pulse.
            mf = 2*sqrt(k)*cos(2*pi*delta*tt) + 1;
            pulses{jj} = pulse(:) .* mf(:) * b1_max_ONres;
     
            
    end
    
end

if nflip==1
    pulses = cell2mat(pulses);
end

% Output modulation function and base pulse if nargout > 1. At the moment
% this only works for the last pulse - use with nflip = 1.
if nargout>1
    pulseBB = pulse(:) * b1_max_ONres;
    pulseMF = mf(:);
    switch nband
        case 2
            pulse_weights = [0 1 2*sqrt(k)];
        case 3
            pulse_weights = [sqrt(k) 1 sqrt(k)];
    end
end

end
