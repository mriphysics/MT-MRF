%%% Simulates the steady-state signal for a bSSFP sequence with multiple
%%% single-band (1B) and multi-band (MB) pulses. Daniel West 2020

function SigDiff = ssSSFP_ihMT_Pulses_2BP(alpha,TR,Dur,dphi,T1x,T2x,M0,K,Delta,TBW,nMB,n1B,B1rms,ReadPoint)

%% Change dimensions of inputs.

TR = TR*1e-3; Dur = Dur*1e-3; Delta = Delta*1e3;

%% Pulse generation.

flips = d2r(alpha);
pulses = {};
[pulses{1},~,~,~,~] = gen_CSMT_pulse_Diffnp(flips,Dur,TR,B1rms,Delta,3,nMB,n1B,'sigma',2);
[pulses{2},~,~,~,~] = gen_CSMT_pulse_Diffnp(flips,Dur,TR,B1rms,Delta,2,nMB,n1B,'sigma',2);

dt = 6.4e-6;
ff = linspace(-20e3,20e3,1000)';
df = ff(2)-ff(1);
nt = length(pulses{1});
tt = dt*(1:nt);
F = exp(-1i*2*pi*ff*tt)*(dt*1e3)/sqrt(numel(ff));

b1sqrd_tau = {}; b1sqrd = {};

band_ix = {};
bw = 2e3;
band_ix{1} = find((ff>-(Delta+bw/2))&(ff<-(Delta-bw/2)));
band_ix{2} = find((ff>-bw/2)&(ff<bw/2));
band_ix{3} = find((ff>(Delta-bw/2))&(ff<(Delta+bw/2)));

for jj=1:2
    
    pwr_spec = abs(F*pulses{jj}).^2;
    
    b1sqrd_tau{jj} = zeros([1 3]);
    for kk=1:3
        b1sqrd_tau{jj}(kk) = sum(pwr_spec(band_ix{kk}))*df;
    end
    tau = 1e3*dt*nt;
    b1sqrd{jj} = b1sqrd_tau{jj}/tau;
    
end

%% Parameters for signal calculation.

B1ms_1B = [0 b1sqrd{1}(1,2) 0]; B1ms_3B = b1sqrd{1}; B1ms_2B = b1sqrd{2};
df_1B = [0 0 0]; df_2B = [0 0 Delta]; df_3B = [-Delta 0 Delta]; % 1B and 2B do not really matter as power will be zero in bands not used.

Gamma = 267.5221;
nd = length(B1ms_1B); % 3 regardless of 1B/2B/3B.

[G_3B,wloc_3B] = SuperLorentzian_LSint(T2x(2),df_3B);
[G_2B,wloc_2B] = SuperLorentzian_LSint(T2x(2),df_2B);
[G_1B,wloc_1B] = SuperLorentzian_LSint(T2x(2),df_1B);
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

T1 = [T1x(1:2) T1x(2) T1x(3)]; % Replicate semisolid Zeeman T1. 
R1 = 1./T1;

% Relaxation/exchange matrices.
E = diag([-R2a_mod -R2a_mod -R1]);
KK = blkdiag(0,0,[[-K*M0s K*M0f K*M0f];[K*(1-delta)*M0s -K*M0f 0];[K*delta*M0s 0 -K*M0f]],0);
O  = diag([-1i*dphi/TR 1i*dphi/TR 0 0 0 0]);
Lambda = E+KK+O; % Matrix exponential of this will form Smat.
C = [0 0 R1(1)*M0f R1(2)*(1-delta)*M0s R1(3)*delta*M0s 0].';

Dphase = diag([-1 -1 1 1 1 1]); % Phase alternation matrix.
% Spoiling matrix would be: Dphase = diag([0 0 1 1 1 1]);

S = [1 0 0 0 0 0]; % Component to return.

%% Perform eigenvector decomposition. Use augmented matrix approach.

LC = cat(1,[Lambda C],zeros(1,7));

AlphaPrime_1B = blkdiag(AMat_1B,1);
AlphaPrime_2B = blkdiag(AMat_2B,1);
AlphaPrime_3B = blkdiag(AMat_3B,1); % Post-expm, so 1 is added rather than 0.

Dphase = blkdiag(Dphase,1);

RDS_1B = AlphaPrime_1B * Dphase * expm(LC*TR);
RDS_2B = AlphaPrime_2B * Dphase * expm(LC*TR);
RDS_3B = AlphaPrime_3B * Dphase * expm(LC*TR);

SRD_1B = expm(LC*TR) * AlphaPrime_1B * Dphase;
SRD_2B = expm(LC*TR) * AlphaPrime_2B * Dphase;
SRD_3B = expm(LC*TR) * AlphaPrime_3B * Dphase;

S_Half = expm(LC*(TR/2));

switch ReadPoint
    
    case 'kernel_end'
        
        X_Total_2B = SRD_1B^n1B * SRD_2B^nMB;
        X_Total_3B = SRD_1B^n1B * SRD_3B^nMB;
        [V_2B,~,~] = eig(X_Total_2B);
        [V_3B,~,~] = eig(X_Total_3B);
        
        MssAll_2B = V_2B(1:end-1,end)/V_2B(end,end);
        MssAll_3B = V_3B(1:end-1,end)/V_3B(end,end); % Normalise to the last component so it stays as 1.
        
        Mss_2B = S*MssAll_2B;
        Mss_3B = S*MssAll_3B;
        
        SigDiff = -abs(Mss_2B - Mss_3B)/sqrt(TR);
        
    case 'after_first_SB'
        
        X_Total_2B = AlphaPrime_1B * Dphase * SRD_2B^nMB * SRD_1B^(n1B-1) * expm(LC*TR);
        X_Total_3B = AlphaPrime_1B * Dphase * SRD_3B^nMB * SRD_1B^(n1B-1) * expm(LC*TR);
        [V_2B,~,~] = eig(X_Total_2B);
        [V_3B,~,~] = eig(X_Total_3B);
        
        MssAll_2B = V_2B(1:end-1,end)/V_2B(end,end);
        MssAll_3B = V_3B(1:end-1,end)/V_3B(end,end);
        
        Mss_2B = S*MssAll_2B;
        Mss_3B = S*MssAll_3B;
        
        SigDiff = -abs(Mss_2B - Mss_3B)/sqrt(TR);
        
    case 'after_last_SB'
        
        X_Total_2B = RDS_1B^n1B * RDS_2B^nMB;
        X_Total_3B = RDS_1B^n1B * RDS_3B^nMB;
        [V_2B,~,~] = eig(X_Total_2B);
        [V_3B,~,~] = eig(X_Total_3B);
        
        MssAll_2B = V_2B(1:end-1,end)/V_2B(end,end);
        MssAll_3B = V_3B(1:end-1,end)/V_3B(end,end);
        
        Mss_2B = S*MssAll_2B;
        Mss_3B = S*MssAll_3B;
        
        SigDiff = -abs(Mss_2B - Mss_3B)/sqrt(TR);
        
    case 'SB_period_middle'
        
        X_Total_2B = S_Half * AlphaPrime_1B * Dphase * SRD_1B^((n1B/2)-1) * SRD_2B^nMB * SRD_1B^(n1B/2) * S_Half;
        X_Total_3B = S_Half * AlphaPrime_1B * Dphase * SRD_1B^((n1B/2)-1) * SRD_3B^nMB * SRD_1B^(n1B/2) * S_Half;
        [V_2B,~,~] = eig(X_Total_2B);
        [V_3B,~,~] = eig(X_Total_3B);
        
        MssAll_2B = V_2B(1:end-1,end)/V_2B(end,end);
        MssAll_3B = V_3B(1:end-1,end)/V_3B(end,end);
        
        Mss_2B = S*MssAll_2B;
        Mss_3B = S*MssAll_3B;
        
        SigDiff = -abs(Mss_2B - Mss_3B)/sqrt(TR);
        
    case 'average_diff'
        
        X_Total_2B_1 = AlphaPrime_1B * Dphase * SRD_2B^nMB * SRD_1B^(n1B-1) * expm(LC*TR);
        X_Total_3B_1 = AlphaPrime_1B * Dphase * SRD_3B^nMB * SRD_1B^(n1B-1) * expm(LC*TR);
        X_Total_2B_2 = S_Half * AlphaPrime_1B * Dphase * SRD_1B^((n1B/2)-1) * SRD_2B^nMB * SRD_1B^(n1B/2) * S_Half;
        X_Total_3B_2 = S_Half * AlphaPrime_1B * Dphase * SRD_1B^((n1B/2)-1) * SRD_3B^nMB * SRD_1B^(n1B/2) * S_Half;
        
        [V_2B_1,~,~] = eig(X_Total_2B_1);
        [V_3B_1,~,~] = eig(X_Total_3B_1);
        [V_2B_2,~,~] = eig(X_Total_2B_2);
        [V_3B_2,~,~] = eig(X_Total_3B_2);
        
        MssAll_2B_1 = V_2B_1(1:end-1,end)/V_2B_1(end,end);
        MssAll_3B_1 = V_3B_1(1:end-1,end)/V_3B_1(end,end);
        MssAll_2B_2 = V_2B_2(1:end-1,end)/V_2B_2(end,end);
        MssAll_3B_2 = V_3B_2(1:end-1,end)/V_3B_2(end,end);
        
        Mss_2B_1 = S*MssAll_2B_1;
        Mss_3B_1 = S*MssAll_3B_1;
        Mss_2B_2 = S*MssAll_2B_2;
        Mss_3B_2 = S*MssAll_3B_2;
        
        SigDiff_1 = -abs(Mss_2B_1 - Mss_3B_1);
        SigDiff_2 = -abs(Mss_2B_2 - Mss_3B_2);
        
        SigDiff = ((SigDiff_1 + SigDiff_2)/2)/sqrt(TR);
        
    case 'MRF'
        
        X_Total_2Bpos = (RDS_1B)*(RDS_2B)^nMB*(RDS_1B)^n1B*(RDS_3B)^nMB*(RDS_1B)^(n1B-1);
        X_Total_3Bpos = (RDS_1B)*(RDS_3B)^nMB*(RDS_1B)^n1B*(RDS_2B)^nMB*(RDS_1B)^(n1B-1);
        
        [V_2Bpos,~,~] = eig(X_Total_2Bpos);
        [V_3Bpos,~,~] = eig(X_Total_3Bpos);
        
        MssAll_2Bpos = V_2Bpos(1:end-1,end)/V_2Bpos(end,end);
        MssAll_3Bpos = V_3Bpos(1:end-1,end)/V_3Bpos(end,end);
       
        Mss_2Bpos = S*MssAll_2Bpos;
        Mss_3Bpos = S*MssAll_3Bpos;
        
        SigDiff = -abs(Mss_2Bpos - Mss_3Bpos)/sqrt(TR);

end

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