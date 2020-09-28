%% Dictionary generation script for fitting. Daniel West 2020.

close all; clear all; clc;

scan_type = 'in_vivo';

%% Generate parameter combinations.

switch scan_type
    % Commented values here are used for 2 bound pool model fit.
    case 'in_vivo'
        n_increments = 100;
        T1f = linspace(0.2,4,n_increments);
        delta = 1; %0.5;
        M0s = linspace(0,0.2,n_increments);
        T1s = 0.2;
        T1d = linspace(0,8e-3,n_increments);
        T2f = 69e-3; %80e-3;
        T2s = 7.5e-6; %12.5e-6;
        k = 50; %65; 
        B0_var = 0; % Assume no off-resonance.
    case 'phantom'
        n_increments = 32;
        T1f = linspace(0.2,4,n_increments);
        delta = 1;
        M0s = linspace(0,0.3,n_increments);
        T1s = 0.3;
        T1d = linspace(0,25e-3,n_increments);
        T2f = 120e-3;
        T2s = 9e-6;
        k = linspace(40,100,n_increments);
        B0_var = 0;
end

T12x_cv = combvec(T1f,T1s,T1d,T2f,T2s).'; T12x_initial_size = size(T12x_cv,1);
% Remove entries where T2a >= T1a (as per Hilbert et al.)
for ii = 1:T12x_initial_size
    if T12x_cv(ii,4) >= T12x_cv(ii,1)
        T12x_cv(ii,:) = NaN;
    end
end
T12x_cv(any(isnan(T12x_cv),2),:) = []; T12x = unique(T12x_cv,'rows'); % Remove repeated rows (if any).

M0_cv = combvec(M0s,delta,k).'; M0_initial_size = size(M0_cv,1);
% If M0b = 0, set f and K to 0.
for ii = 1:M0_initial_size
   if M0_cv(ii,1) == 0
        M0_cv(ii,2:3) = 0;
   end
end

M0f = 1-M0_cv(:,1); M0_comp = unique([M0f, M0_cv],'rows'); % M0f is first dim, M0s is second and delta is third.

noMT_idx = find(M0_comp(:,2) == 0); M0_noMT = M0_comp(noMT_idx,:);
MT_idx = find(M0_comp(:,3) == 0 & M0_comp(:,2) ~= 0); M0_MT = M0_comp(MT_idx,:);
ihMT_idx = find(M0_comp(:,3) ~= 0); M0_ihMT = M0_comp(ihMT_idx,:);

% Combine all input parameters into a single matrix.
Var_Params_noMT = unique(combvec(M0_noMT.',T12x.',B0_var).','rows');
Var_Params_MT = unique(combvec(M0_MT.',T12x.',B0_var).','rows');
Var_Params_ihMT = unique(combvec(M0_ihMT.',T12x.',B0_var).','rows');

%% Set-up sequence.

flips = d2r(29.51); TR = 5.33e-3; Dur = 2.51e-3; Delta = 8058.48; n1B = 300; nMB = 300;
TBW = 2; B1rms = 4;
  
% Generate pulses.
pulses = {};
[pulses{1,1},~,~,~,~] = gen_CSMT_pulse_Diffnp(flips,Dur,TR,B1rms,Delta,2,nMB,n1B,'sigma',2);
[pulses{1,2},~,~,~,~] = gen_CSMT_pulse_Diffnp(flips,Dur,TR,B1rms,Delta,3,nMB,n1B,'sigma',2);
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
for jj = 1:2
    pwr_spec = abs(F*pulses{1,jj}).^2;
    b1sqrd_tau{1,jj} = zeros([1 3]);
    for kk=1:3
        b1sqrd_tau{1,jj}(kk) = sum(pwr_spec(band_ix{kk}))*df;
    end
    tau = 1e3*dt*nt;
    b1sqrd{1,jj} = b1sqrd_tau{1,jj}/tau;
end

np_Total = (2*nMB)+(2*n1B); % Total number of pulses.

for ii = 1:length(T2s)
    % Pre-compute lineshape for different T2b - removed loop because T2b fixed.
    df_1B = [0 0 0]; df_2B = [0 0 Delta]; df_3B = [-Delta 0 Delta];
    [G_3B(ii,:),wloc_3B(ii)] = SuperLorentzian_LSint(T2s(ii),df_3B);
    [G_2B(ii,:),wloc_2B(ii)] = SuperLorentzian_LSint(T2s(ii),df_2B);
    [G_1B(ii,:),wloc_1B(ii)] = SuperLorentzian_LSint(T2s(ii),df_1B);
end

% Calculates signals for one sequence "period".
Mss_Total_noMT = zeros(size(Var_Params_noMT,1),np_Total);
Mss_Total_MT = zeros(size(Var_Params_MT,1),np_Total);
Mss_Total_ihMT = zeros(size(Var_Params_ihMT,1),np_Total);

Dict_size = GetSize(Mss_Total_noMT) + GetSize(Mss_Total_MT) + GetSize(Mss_Total_ihMT);
disp(['Total dictionary size is: ',num2str(Dict_size/1e9),'GB.'])

%% Signal generation.

delete(gcp('nocreate')); c = parcluster('local'); c.NumWorkers = 16; parpool(c, c.NumWorkers);

% Form M0 matrix for signal generation function.
M0_noMT_fit = [Var_Params_noMT(:,1) , Var_Params_noMT(:,2).*(1-Var_Params_noMT(:,3)) , Var_Params_noMT(:,2).*Var_Params_noMT(:,3)]; 
M0_MT_fit = [Var_Params_MT(:,1) , Var_Params_MT(:,2).*(1-Var_Params_MT(:,3)) , Var_Params_MT(:,2).*Var_Params_MT(:,3)]; 
M0_ihMT_fit = [Var_Params_ihMT(:,1) , Var_Params_ihMT(:,2).*(1-Var_Params_ihMT(:,3)) , Var_Params_ihMT(:,2).*Var_Params_ihMT(:,3)]; 

% No MT dictionary creation.
tic
parfor mm = 1:size(Var_Params_noMT,1)
    idx = find(T2s(:) == Var_Params_noMT(mm,9));
    Mss_Total_noMT(mm,:) = Dictionary_function_CSS(flips,TR,Dur,Var_Params_noMT(mm,10),Var_Params_noMT(mm,5:7),Var_Params_noMT(mm,8:9),M0_noMT_fit(mm,:),Var_Params_noMT(mm,4),Delta,TBW,nMB,n1B,b1sqrd,G_1B(idx,:),G_2B(idx,:),G_3B(idx,:),wloc_1B(idx),wloc_2B(idx),wloc_3B(idx));
end
toc

% MT dictionary creation.
tic
parfor mm = 1:size(Var_Params_MT,1)
    idx = find(T2s(:) == Var_Params_MT(mm,9));
    Mss_Total_MT(mm,:) = Dictionary_function_CSS(flips,TR,Dur,Var_Params_MT(mm,10),Var_Params_MT(mm,5:7),Var_Params_MT(mm,8:9),M0_MT_fit(mm,:),Var_Params_MT(mm,4),Delta,TBW,nMB,n1B,b1sqrd,G_1B(idx,:),G_2B(idx,:),G_3B(idx,:),wloc_1B(idx),wloc_2B(idx),wloc_3B(idx));
end
toc

% ihMT dictionary creation.
tic
parfor mm = 1:size(Var_Params_ihMT,1)
    idx = find(T2s(:) == Var_Params_ihMT(mm,9));
    Mss_Total_ihMT(mm,:) = Dictionary_function_CSS(flips,TR,Dur,Var_Params_ihMT(mm,10),Var_Params_ihMT(mm,5:7),Var_Params_ihMT(mm,8:9),M0_ihMT_fit(mm,:),Var_Params_ihMT(mm,4),Delta,TBW,nMB,n1B,b1sqrd,G_1B(idx,:),G_2B(idx,:),G_3B(idx,:),wloc_1B(idx),wloc_2B(idx),wloc_3B(idx));
end
toc

Mss_AllSigs = [Mss_Total_noMT ; Mss_Total_MT ; Mss_Total_ihMT];
VP_AllSigs = [Var_Params_noMT ; Var_Params_MT ; Var_Params_ihMT];

delete(gcp('nocreate')); clearvars -except Mss_AllSigs VP_AllSigs
