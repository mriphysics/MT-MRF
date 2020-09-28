%%% Simulates the steady-state signal for a bSSFP sequence with multiple
%%% single-band (1B) and multi-band (MB) pulses. Daniel West 2020

function Mss_Total = Dictionary_function_CSS(flips,TR,Dur,dphi,T1x,T2x,M0,K,Delta,TBW,nMB,n1B,b1sqrd,G_1B,G_2B,G_3B,wloc_1B,wloc_2B,wloc_3B)

%% Parameters for signal calculation.

B1ms_1B = [0 b1sqrd{1}(1,2) 0]; B1ms_2B = b1sqrd{1}; B1ms_3B = b1sqrd{2};
df_1B = [0 0 0]; df_2B = [0 0 Delta]; df_3B = [-Delta 0 Delta];

Gamma = 267.52219; nd = length(B1ms_1B);

W_3B = pi*Gamma^2*B1ms_3B.*G_3B;
W_2B = pi*Gamma^2*B1ms_2B.*G_2B; 
W_1B = pi*Gamma^2*B1ms_1B.*G_1B;

%% Evolution during RF - flip and saturate. Two separate terms needed for 1B and MB.

WPulse_1B = zeros(3,3); WPulse_2B = zeros(3,3); WPulse_3B = zeros(3,3);
for ii=1:nd
    D = 2*pi*df_1B(ii)/wloc_1B;
    WPulse_1B = WPulse_1B + W_1B(ii)*([[-1 0 0];[0 -1 D];[0 D -D^2]]);
end
for ii=1:nd
    D = 2*pi*df_2B(ii)/wloc_2B;
    WPulse_2B = WPulse_2B + W_2B(ii)*([[-1 0 0];[0 -1 D];[0 D -D^2]]);
end
for ii=1:nd
    D = 2*pi*df_3B(ii)/wloc_3B;
    WPulse_3B = WPulse_3B + W_3B(ii)*([[-1 0 0];[0 -1 D];[0 D -D^2]]);
end

% Assume RF is rotation around x-axis. Phase doesn't matter here.
% Pick the sum(alpha) assuming central is the only non-zero frequency band.
T = RF_Rot(sum(flips),0);
AMat_1B = blkdiag(T,expm(WPulse_1B*Dur));
AMat_2B = blkdiag(T,expm(WPulse_2B*Dur));
AMat_3B = blkdiag(T,expm(WPulse_3B*Dur));

%% Evolution during free precession. 

% Modify T2 as per Bieri et al. 
R1a = 1/T1x(1); R2a = 1/T2x(1); 
T_RFE = 1.20*(Dur/TBW);
Eta = 0.68-0.125*(1+(T_RFE/TR))*(R1a/R2a);
R2a_mod = (1-Eta*(T_RFE/TR))*R2a;

M0f = M0(1);
M0s = M0(2)+M0(3); % Total semisolid fraction.
delta = M0(3)/M0s; % Fraction of semisolid pool that is dipolar coupled.

% Reset f value from NaN.
if M0s == 0; delta = 0; end

T1 = [T1x(1:2) T1x(2) T1x(3)]; % Replicate semisolid Zeeman T1. 
R1 = 1./T1;

% Reset R1D value from Inf.
if R1(4) == Inf; R1(4) = 1e12; end

% Relaxation/exchange matrices.
E = diag([-R2a_mod -R2a_mod -R1]);
KK = blkdiag(0,0,[[-K*M0s K*M0f K*M0f];[K*(1-delta)*M0s -K*M0f 0];[K*delta*M0s 0 -K*M0f]],0);
O  = diag([-1i*dphi/TR 1i*dphi/TR 0 0 0 0]);
Lambda = E+KK+O; % Matrix exponential of this will form Smat.
C = [0 0 R1(1)*M0f R1(2)*(1-delta)*M0s R1(3)*delta*M0s 0].';

Dphase = diag([-1 -1 1 1 1 1]); % Phase alternation matrix.
% Spoiling matrix would be: Dphase = diag([0 0 1 1 1 1]);

%% Perform eigenvector decomposition. Use augmented matrix approach.

LC = cat(1,[Lambda C],zeros(1,7));

AlphaPrime_1B = blkdiag(AMat_1B,1); 
AlphaPrime_2B = blkdiag(AMat_2B,1); 
AlphaPrime_3B = blkdiag(AMat_3B,1); % Post-expm, so 1 is added rather than 0.

Dphase = blkdiag(Dphase,1);

RDS_1B = AlphaPrime_1B * Dphase * expm(LC*TR); 
RDS_2B = AlphaPrime_2B * Dphase * expm(LC*TR);
RDS_3B = AlphaPrime_3B * Dphase * expm(LC*TR);

%%% Signal calculation. %%%

% First reference point.
X = RDS_2B * (RDS_1B)^n1B * (RDS_3B)^nMB * (RDS_1B)^n1B * (RDS_2B)^(nMB-1);
[V,~,~] = eig(X);

% Create array to hold all state vectors.
Mstates = zeros([size(X,1) 2*(n1B+nMB)]);

Mstates(:,1) = V(1:end,end)/V(end,end);

Counter = 2;
% 2-band period.
for ii = 2:nMB
   Mstates(:,Counter) = RDS_2B * Mstates(:,Counter-1);
   Counter=Counter+1;
end
% 1-band period.
for ii = 1:n1B
   Mstates(:,Counter) = RDS_1B * Mstates(:,Counter-1);
   Counter=Counter+1;
end
% 3-band period.
for ii = 1:nMB
   Mstates(:,Counter) = RDS_3B * Mstates(:,Counter-1);
   Counter=Counter+1;
end
% 1-band period.
for ii = 1:n1B
   Mstates(:,Counter) = RDS_1B * Mstates(:,Counter-1);
   Counter=Counter+1;
end

Mss_Total = Mstates(1,:);

%%% Deals with rotation. %%%
    function Tap = RF_Rot(a,p)
        Tap = zeros([3 3]);
        Tap(1) = cos(a/2).^2;
        Tap(2) = exp(-2*1i*p)*(sin(a/2)).^2;
        Tap(3) = -0.5*1i*exp(-1i*p)*sin(a);
        Tap(4) = conj(Tap(2));
        Tap(5) = Tap(1);
        Tap(6) = 0.5*1i*exp(1i*p)*sin(a);
        Tap(7) = -1i*exp(1i*p)*sin(a);
        Tap(8) = 1i*exp(-1i*p)*sin(a);
        Tap(9) = cos(a);
    end

end