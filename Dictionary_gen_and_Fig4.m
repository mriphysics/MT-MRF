%% MTF dictionary generation script. Daniel West 2020.

close all; clear all; clc;

%% Generate different parameter combinations.

% Specify ranges for tissue parameters.
T1f = 0.2:0.5:4; T1s = 0.1:0.2:1; T1d = 2e-3:4e-3:25e-3;
T2f = 50e-3:70e-3:500e-3; T2s = 5e-6:3e-6:20e-6;
M0s = 0:0.05:0.25; delta = 0:0.2:1; k = 40:20:100;
B0_var = -pi:pi/3:pi;

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

% Form M0 matrix and remove repeated rows (if any).
M0_noMT = unique([M0_noMT(:,1) , M0_noMT(:,2).*(1-M0_noMT(:,3)) , M0_noMT(:,2).*M0_noMT(:,3)],'rows'); 
M0_MT = unique([M0_MT(:,1) , M0_MT(:,2).*(1-M0_MT(:,3)) , M0_MT(:,2).*M0_MT(:,3)],'rows'); 
M0_ihMT = unique([M0_ihMT(:,1) , M0_ihMT(:,2).*(1-M0_ihMT(:,3)) , M0_ihMT(:,2).*M0_ihMT(:,3)],'rows'); 

% Combine all input parameters into a single matrix.
Var_Params_noMT = unique(combvec(M0_noMT.',T12x.',0,B0_var).','rows');
Var_Params_MT = unique(combvec(M0_MT.',T12x.',k,B0_var).','rows');
Var_Params_ihMT = unique(combvec(M0_ihMT.',T12x.',k,B0_var).','rows');

% Remove infeasible elements for each tissue.
for ii = 1:size(Var_Params_noMT,1)
    if Var_Params_noMT(ii,4) < 2.5 || Var_Params_noMT(ii,7) < 300e-3
        Var_Params_noMT(ii,:) = NaN;
    end
end
Var_Params_noMT(any(isnan(Var_Params_noMT),2),:) = [];
for ii = 1:size(Var_Params_MT,1)
    if Var_Params_MT(ii,4) > 2.5 || Var_Params_MT(ii,7) > 300e-3 
        Var_Params_MT(ii,:) = NaN;
    end
end
Var_Params_MT(any(isnan(Var_Params_MT),2),:) = [];
for ii = 1:size(Var_Params_ihMT,1)
    if Var_Params_ihMT(ii,4) > 1.5 || Var_Params_ihMT(ii,7) > 150e-3
        Var_Params_ihMT(ii,:) = NaN;
    end
end
Var_Params_ihMT(any(isnan(Var_Params_ihMT),2),:) = [];

%% Specify sequence details.

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

% Pre-compute lineshape for different T2b - removed loop because T2b fixed.
df_1B = [0 0 0]; df_2B = [0 0 Delta]; df_3B = [-Delta 0 Delta];
for ii = 1:length(T2s)
    [G_3B(ii,:),wloc_3B(ii)] = SuperLorentzian_LSint(T2s(ii),df_3B);
    [G_2B(ii,:),wloc_2B(ii)] = SuperLorentzian_LSint(T2s(ii),df_2B);
    [G_1B(ii,:),wloc_1B(ii)] = SuperLorentzian_LSint(T2s(ii),df_1B);
end

% Calculates signals for one sequence "period".
Mss_Total_noMT = zeros(size(Var_Params_noMT,1),np_Total); whos Mss_Total_noMT;
Mss_Total_MT = zeros(size(Var_Params_MT,1),np_Total); whos Mss_Total_MT;
Mss_Total_ihMT = zeros(size(Var_Params_ihMT,1),np_Total); whos Mss_Total_ihMT;

%% Simulate signals for different cparameter combinations.

delete(gcp('nocreate')); c = parcluster('local'); c.NumWorkers = 16; parpool(c, c.NumWorkers);

% No MT dictionary creation.
tic
parfor mm = 1:size(Var_Params_noMT,1)
    if Var_Params_noMT(mm,8) == T2s(1)
        Mss_Total_noMT(mm,:) = Dictionary_function_CSS(flips,TR,Dur,Var_Params_noMT(mm,10),Var_Params_noMT(mm,4:6),Var_Params_noMT(mm,7:8),Var_Params_noMT(mm,1:3),Var_Params_noMT(mm,9),Delta,TBW,nMB,n1B,b1sqrd,G_1B(1,:),G_2B(1,:),G_3B(1,:),wloc_1B(1),wloc_2B(1),wloc_3B(1));
    end
    if Var_Params_noMT(mm,8) == T2s(2)
        Mss_Total_noMT(mm,:) = Dictionary_function_CSS(flips,TR,Dur,Var_Params_noMT(mm,10),Var_Params_noMT(mm,4:6),Var_Params_noMT(mm,7:8),Var_Params_noMT(mm,1:3),Var_Params_noMT(mm,9),Delta,TBW,nMB,n1B,b1sqrd,G_1B(2,:),G_2B(2,:),G_3B(2,:),wloc_1B(2),wloc_2B(2),wloc_3B(2));
    end
    if Var_Params_noMT(mm,8) == T2s(3)
        Mss_Total_noMT(mm,:) = Dictionary_function_CSS(flips,TR,Dur,Var_Params_noMT(mm,10),Var_Params_noMT(mm,4:6),Var_Params_noMT(mm,7:8),Var_Params_noMT(mm,1:3),Var_Params_noMT(mm,9),Delta,TBW,nMB,n1B,b1sqrd,G_1B(3,:),G_2B(3,:),G_3B(3,:),wloc_1B(3),wloc_2B(3),wloc_3B(3));
    end
    if Var_Params_noMT(mm,8) == T2s(4)
        Mss_Total_noMT(mm,:) = Dictionary_function_CSS(flips,TR,Dur,Var_Params_noMT(mm,10),Var_Params_noMT(mm,4:6),Var_Params_noMT(mm,7:8),Var_Params_noMT(mm,1:3),Var_Params_noMT(mm,9),Delta,TBW,nMB,n1B,b1sqrd,G_1B(4,:),G_2B(4,:),G_3B(4,:),wloc_1B(4),wloc_2B(4),wloc_3B(4));
    end
    if Var_Params_noMT(mm,8) == T2s(5)
        Mss_Total_noMT(mm,:) = Dictionary_function_CSS(flips,TR,Dur,Var_Params_noMT(mm,10),Var_Params_noMT(mm,4:6),Var_Params_noMT(mm,7:8),Var_Params_noMT(mm,1:3),Var_Params_noMT(mm,9),Delta,TBW,nMB,n1B,b1sqrd,G_1B(5,:),G_2B(5,:),G_3B(5,:),wloc_1B(5),wloc_2B(5),wloc_3B(5));
    end
    if Var_Params_noMT(mm,8) == T2s(6)
        Mss_Total_noMT(mm,:) = Dictionary_function_CSS(flips,TR,Dur,Var_Params_noMT(mm,10),Var_Params_noMT(mm,4:6),Var_Params_noMT(mm,7:8),Var_Params_noMT(mm,1:3),Var_Params_noMT(mm,9),Delta,TBW,nMB,n1B,b1sqrd,G_1B(6,:),G_2B(6,:),G_3B(6,:),wloc_1B(6),wloc_2B(6),wloc_3B(6));
    end
end
toc

% MT dictionary creation.
tic
parfor mm = 1:size(Var_Params_MT,1)
    if Var_Params_MT(mm,8) == T2s(1)
        Mss_Total_MT(mm,:) = Dictionary_function_CSS(flips,TR,Dur,Var_Params_MT(mm,10),Var_Params_MT(mm,4:6),Var_Params_MT(mm,7:8),Var_Params_MT(mm,1:3),Var_Params_MT(mm,9),Delta,TBW,nMB,n1B,b1sqrd,G_1B(1,:),G_2B(1,:),G_3B(1,:),wloc_1B(1),wloc_2B(1),wloc_3B(1));
    end
    if Var_Params_MT(mm,8) == T2s(2)
        Mss_Total_MT(mm,:) = Dictionary_function_CSS(flips,TR,Dur,Var_Params_MT(mm,10),Var_Params_MT(mm,4:6),Var_Params_MT(mm,7:8),Var_Params_MT(mm,1:3),Var_Params_MT(mm,9),Delta,TBW,nMB,n1B,b1sqrd,G_1B(2,:),G_2B(2,:),G_3B(2,:),wloc_1B(2),wloc_2B(2),wloc_3B(2));
    end
    if Var_Params_MT(mm,8) == T2s(3)
        Mss_Total_MT(mm,:) = Dictionary_function_CSS(flips,TR,Dur,Var_Params_MT(mm,10),Var_Params_MT(mm,4:6),Var_Params_MT(mm,7:8),Var_Params_MT(mm,1:3),Var_Params_MT(mm,9),Delta,TBW,nMB,n1B,b1sqrd,G_1B(3,:),G_2B(3,:),G_3B(3,:),wloc_1B(3),wloc_2B(3),wloc_3B(3));
    end
    if Var_Params_MT(mm,8) == T2s(4)
        Mss_Total_MT(mm,:) = Dictionary_function_CSS(flips,TR,Dur,Var_Params_MT(mm,10),Var_Params_MT(mm,4:6),Var_Params_MT(mm,7:8),Var_Params_MT(mm,1:3),Var_Params_MT(mm,9),Delta,TBW,nMB,n1B,b1sqrd,G_1B(4,:),G_2B(4,:),G_3B(4,:),wloc_1B(4),wloc_2B(4),wloc_3B(4));
    end
    if Var_Params_MT(mm,8) == T2s(5)
        Mss_Total_MT(mm,:) = Dictionary_function_CSS(flips,TR,Dur,Var_Params_MT(mm,10),Var_Params_MT(mm,4:6),Var_Params_MT(mm,7:8),Var_Params_MT(mm,1:3),Var_Params_MT(mm,9),Delta,TBW,nMB,n1B,b1sqrd,G_1B(5,:),G_2B(5,:),G_3B(5,:),wloc_1B(5),wloc_2B(5),wloc_3B(5));
    end
    if Var_Params_MT(mm,8) == T2s(6)
        Mss_Total_MT(mm,:) = Dictionary_function_CSS(flips,TR,Dur,Var_Params_MT(mm,10),Var_Params_MT(mm,4:6),Var_Params_MT(mm,7:8),Var_Params_MT(mm,1:3),Var_Params_MT(mm,9),Delta,TBW,nMB,n1B,b1sqrd,G_1B(6,:),G_2B(6,:),G_3B(6,:),wloc_1B(6),wloc_2B(6),wloc_3B(6));
    end
end
toc

% ihMT dictionary creation.
tic
parfor mm = 1:size(Var_Params_ihMT,1)
    if Var_Params_ihMT(mm,8) == T2s(1)
        Mss_Total_ihMT(mm,:) = Dictionary_function_CSS(flips,TR,Dur,Var_Params_ihMT(mm,10),Var_Params_ihMT(mm,4:6),Var_Params_ihMT(mm,7:8),Var_Params_ihMT(mm,1:3),Var_Params_ihMT(mm,9),Delta,TBW,nMB,n1B,b1sqrd,G_1B(1,:),G_2B(1,:),G_3B(1,:),wloc_1B(1),wloc_2B(1),wloc_3B(1));
    end
    if Var_Params_ihMT(mm,8) == T2s(2)
        Mss_Total_ihMT(mm,:) = Dictionary_function_CSS(flips,TR,Dur,Var_Params_ihMT(mm,10),Var_Params_ihMT(mm,4:6),Var_Params_ihMT(mm,7:8),Var_Params_ihMT(mm,1:3),Var_Params_ihMT(mm,9),Delta,TBW,nMB,n1B,b1sqrd,G_1B(2,:),G_2B(2,:),G_3B(2,:),wloc_1B(2),wloc_2B(2),wloc_3B(2));
    end
    if Var_Params_ihMT(mm,8) == T2s(3)
        Mss_Total_ihMT(mm,:) = Dictionary_function_CSS(flips,TR,Dur,Var_Params_ihMT(mm,10),Var_Params_ihMT(mm,4:6),Var_Params_ihMT(mm,7:8),Var_Params_ihMT(mm,1:3),Var_Params_ihMT(mm,9),Delta,TBW,nMB,n1B,b1sqrd,G_1B(3,:),G_2B(3,:),G_3B(3,:),wloc_1B(3),wloc_2B(3),wloc_3B(3));
    end
    if Var_Params_ihMT(mm,8) == T2s(4)
        Mss_Total_ihMT(mm,:) = Dictionary_function_CSS(flips,TR,Dur,Var_Params_ihMT(mm,10),Var_Params_ihMT(mm,4:6),Var_Params_ihMT(mm,7:8),Var_Params_ihMT(mm,1:3),Var_Params_ihMT(mm,9),Delta,TBW,nMB,n1B,b1sqrd,G_1B(4,:),G_2B(4,:),G_3B(4,:),wloc_1B(4),wloc_2B(4),wloc_3B(4));
    end
    if Var_Params_ihMT(mm,8) == T2s(5)
        Mss_Total_ihMT(mm,:) = Dictionary_function_CSS(flips,TR,Dur,Var_Params_ihMT(mm,10),Var_Params_ihMT(mm,4:6),Var_Params_ihMT(mm,7:8),Var_Params_ihMT(mm,1:3),Var_Params_ihMT(mm,9),Delta,TBW,nMB,n1B,b1sqrd,G_1B(5,:),G_2B(5,:),G_3B(5,:),wloc_1B(5),wloc_2B(5),wloc_3B(5));
    end
    if Var_Params_ihMT(mm,8) == T2s(6)
        Mss_Total_ihMT(mm,:) = Dictionary_function_CSS(flips,TR,Dur,Var_Params_ihMT(mm,10),Var_Params_ihMT(mm,4:6),Var_Params_ihMT(mm,7:8),Var_Params_ihMT(mm,1:3),Var_Params_ihMT(mm,9),Delta,TBW,nMB,n1B,b1sqrd,G_1B(6,:),G_2B(6,:),G_3B(6,:),wloc_1B(6),wloc_2B(6),wloc_3B(6));
    end
end
toc

%% Select representative samples from dictionary.

rng('default');
% Truncate dictionary by a certain factor/make dictionary representative.
k_noMT = randperm(size(Var_Params_noMT,1),10000).';
k_MT = randperm(size(Var_Params_MT,1),200000).';
k_ihMT = randperm(size(Var_Params_ihMT,1),200000).';

Tr_Mss_Total_noMT = Mss_Total_noMT(k_noMT,:);
Tr_Mss_Total_MT = Mss_Total_MT(k_MT,:);
Tr_Mss_Total_ihMT = Mss_Total_ihMT(k_ihMT,:);

Tr_Var_Params_noMT = Var_Params_noMT(k_noMT,:);
Tr_Var_Params_MT = Var_Params_MT(k_MT,:);
Tr_Var_Params_ihMT = Var_Params_ihMT(k_ihMT,:);

clear Mss_Total_noMT Mss_Total_MT Mss_Total_ihMT % Saves memory.

% Create truncated dictionary.
Tr_Mss_AllSigs = [repmat(Tr_Mss_Total_noMT,20,1) ; Tr_Mss_Total_MT ; Tr_Mss_Total_ihMT];
Tr_VP_AllSigs = [repmat(Tr_Var_Params_noMT,20,1) ; Tr_Var_Params_MT ; Tr_Var_Params_ihMT];

%% Take SVD of dictionary.

[U,S,V] = svd(Tr_Mss_AllSigs,'econ');
disp(['5 largest singular values: ', num2str(S(1,1)), ', ', num2str(S(2,2)), ', ', num2str(S(3,3)), ', ', num2str(S(4,4)), ', ', num2str(S(5,5)), '.'])

%% Dictionary analysis for different simulated tissues.

% WM parameters.
T1x_WM = [650 1000 6.5]*1e-3; 
T2x_WM = [80e-3, 12.5e-6]; 
f_WM = 0.65; 
K_WM = 65; 
M0b_Varma_WM = (7.3*(1/T1x_WM(1)))/K_WM; 
M0b_WM = M0b_Varma_WM/(1 + M0b_Varma_WM); 
M0f_WM = 1-M0b_WM; 
M0_WM = [M0f_WM M0b_WM*(1-f_WM) M0b_WM*f_WM];
VP_WM = [M0_WM, T1x_WM, T2x_WM, K_WM];
[GWM_3B,wlocWM_3B] = SuperLorentzian_LSint(T2x_WM(2),df_3B);
[GWM_2B,wlocWM_2B] = SuperLorentzian_LSint(T2x_WM(2),df_2B);
[GWM_1B,wlocWM_1B] = SuperLorentzian_LSint(T2x_WM(2),df_1B);

% GM parameters.
T1x_GM = [1100 1000 5.5]*1e-3; 
T2x_GM = [115e-3, 10e-6]; 
f_GM = 0.20; 
K_GM = 50; 
M0b_Varma_GM = (7.3*(1/T1x_GM(1)))/K_GM; 
M0b_GM = M0b_Varma_GM/(1 + M0b_Varma_GM); 
M0f_GM = 1-M0b_GM; 
M0_GM = [M0f_GM M0b_GM*(1-f_GM) M0b_GM*f_GM];
VP_GM = [M0_GM, T1x_GM, T2x_GM, K_GM];
[GGM_3B,wlocGM_3B] = SuperLorentzian_LSint(T2x_GM(2),df_3B);
[GGM_2B,wlocGM_2B] = SuperLorentzian_LSint(T2x_GM(2),df_2B);
[GGM_1B,wlocGM_1B] = SuperLorentzian_LSint(T2x_GM(2),df_1B);

% CSF parmeters.
T1x_CSF = [3000 1000 1000]*1e-3; 
T2x_CSF = [2, 1]; 
f_CSF = 0; 
K_CSF = 0; 
M0b_Varma_CSF = 0;
M0b_CSF = M0b_Varma_CSF/(1 + M0b_Varma_CSF); 
M0f_CSF = 1-M0b_CSF; 
M0_CSF = [M0f_CSF M0b_CSF*(1-f_CSF) M0b_CSF*f_CSF];
VP_CSF = [M0_CSF, T1x_CSF, T2x_CSF, K_CSF];
[GCSF_3B,wlocCSF_3B] = SuperLorentzian_LSint(T2x_CSF(2),df_3B);
[GCSF_2B,wlocCSF_2B] = SuperLorentzian_LSint(T2x_CSF(2),df_2B);
[GCSF_1B,wlocCSF_1B] = SuperLorentzian_LSint(T2x_CSF(2),df_1B);

Rank_no = 1:20; dphi = 0;
x_WM = zeros(length(Rank_no),1200);xprime_WM = zeros(length(Rank_no),1200);
x_GM = zeros(length(Rank_no),1200);xprime_GM = zeros(length(Rank_no),1200);
x_CSF = zeros(length(Rank_no),1200);xprime_CSF = zeros(length(Rank_no),1200);
for ii = 1:length(Rank_no)

    % Generate signals.
    x_WM(ii,:) = (Dictionary_function_CSS(flips,TR,Dur,dphi,VP_WM(4:6),VP_WM(7:8),VP_WM(1:3),VP_WM(9),Delta,TBW,nMB,n1B,b1sqrd,GWM_1B,GWM_2B,GWM_3B,wlocWM_1B,wlocWM_2B,wlocWM_3B));
    w_WM = pinv(V(:,1:Rank_no(ii)))*x_WM(ii,:).';
    xprime_WM(ii,:) = V(:,1:Rank_no(ii))*w_WM;
    x_GM(ii,:) = (Dictionary_function_CSS(flips,TR,Dur,dphi,VP_GM(4:6),VP_GM(7:8),VP_GM(1:3),VP_GM(9),Delta,TBW,nMB,n1B,b1sqrd,GGM_1B,GGM_2B,GGM_3B,wlocGM_1B,wlocGM_2B,wlocGM_3B));
    w_GM = pinv(V(:,1:Rank_no(ii)))*x_GM(ii,:).';
    xprime_GM(ii,:) = V(:,1:Rank_no(ii))*w_GM;
    x_CSF(ii,:) = (Dictionary_function_CSS(flips,TR,Dur,dphi,VP_CSF(4:6),VP_CSF(7:8),VP_CSF(1:3),VP_CSF(9),Delta,TBW,nMB,n1B,b1sqrd,GCSF_1B,GCSF_2B,GCSF_3B,wlocCSF_1B,wlocCSF_2B,wlocCSF_3B));
    w_CSF = pinv(V(:,1:Rank_no(ii)))*x_CSF(ii,:).';
    xprime_CSF(ii,:) = V(:,1:Rank_no(ii))*w_CSF;

end

%% Plot results of analysis for Figure 4.

cm = lines(10);
figure(1); 
labels = {'R = 1','R = 2','R = 3','R = 4','R = 5','R = 6','R = 7','R = 8','R = 9','R = 10'};

a_plot(:,2) = V(:,1); a_plot(:,1) = zeros(1200,1);
b_plot(:,2) = V(:,2); b_plot(:,1) = zeros(1200,1);
c_plot(:,2) = V(:,3); c_plot(:,1) = zeros(1200,1);
d_plot(:,2) = V(:,4); d_plot(:,1) = zeros(1200,1);
e_plot(:,2) = V(:,5); e_plot(:,1) = zeros(1200,1);
f_plot(:,2) = V(:,6); f_plot(:,1) = zeros(1200,1);
g_plot(:,2) = V(:,7); g_plot(:,1) = zeros(1200,1);
h_plot(:,2) = V(:,8); h_plot(:,1) = zeros(1200,1);
i_plot(:,2) = V(:,9); i_plot(:,1) = zeros(1200,1);
j_plot(:,2) = V(:,10); j_plot(:,1) = zeros(1200,1);

t = table(a_plot,b_plot,c_plot,d_plot,e_plot,f_plot,g_plot,h_plot,i_plot,j_plot);

subplot(10,2,1:2:19)
s = stackedplot(t); set(gca,'FontSize',18); s.GridVisible = 'on'; s.XLabel = {'Pulse No.'};
for ii = 1:10
    s.AxesProperties(ii).YLimits = [-0.2 0.2];
    s.AxesProperties(ii).LegendVisible = 'off';
    s.DisplayLabels(ii) = labels(ii);
    s.LineProperties(ii).Color = [[0.6 0.6 0.6];cm(ii,:)];
    s.LineProperties(ii).LineStyle = {'--','-'};
    s.LineProperties(ii).LineWidth = 3;
end
annotation('textbox','String','(a)','FontSize',24,'LineStyle','none','Position',[0.0644,0.905,0.05,0.07])


subplot(10,2,2:2:10)
plot(log10(diag(S)),'-o','LineWidth',3,'MarkerSize',10,'Color','k'); xlim([0 20])
set(gca,'FontSize',16); grid on; grid minor; xlabel('R','FontSize',20); ylabel('Log of Singular Values');
annotation('textbox','String','(b)','FontSize',24,'LineStyle','none','Position',[0.5221875,0.905,0.05,0.069999999999999])

% Compute mean and max residuals.
Maxres_WM = max(((abs(x_WM(5,:)-xprime_WM(5,:)))./abs(x_WM(5,:)))*100);
Maxres_GM = max(((abs(x_GM(5,:)-xprime_GM(5,:)))./abs(x_GM(5,:)))*100);
Maxres_CSF = max(((abs(x_CSF(5,:)-xprime_CSF(5,:)))./abs(x_CSF(5,:)))*100);

Meanres_WM = zeros(1,10); Meanres_GM = zeros(1,10); Meanres_CSF = zeros(1,10);
for ii = 1:10
    Meanres_WM(ii) = mean(((abs(x_WM(ii,:)-xprime_WM(ii,:)))./abs(x_WM(ii,:)))*100);
    Meanres_GM(ii) = mean(((abs(x_GM(ii,:)-xprime_GM(ii,:)))./abs(x_GM(ii,:)))*100);
    Meanres_CSF(ii) = mean(((abs(x_CSF(ii,:)-xprime_CSF(ii,:)))./abs(x_CSF(ii,:)))*100);
end

cm2 = gray(20);
subplot(10,2,14:2:20)
plot(Meanres_WM,'-o','LineWidth',3,'MarkerSize',10,'Color',cm2(1,:)); hold on
plot(Meanres_GM,'-o','LineWidth',3,'MarkerSize',10,'Color',cm2(8,:));
plot(Meanres_CSF,'-o','LineWidth',3,'MarkerSize',10,'Color',cm2(15,:));
set(gca,'FontSize',16); xlabel('R','FontSize',20); ylabel('Mean Residuals (%)');
grid on; grid minor; ll = legend('WM','GM','CSF'); ll.FontSize = 20; ll.Orientation = 'horizontal'; legend boxoff; xlim([1 10])
annotation('textbox','String','(c)','FontSize',24,'LineStyle','none','Position',[0.5221875,0.418,0.05,0.07])