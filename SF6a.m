%%% R2: Numerical simulation study. Part 1 %%%

close all; clear all; clc

%% Generate ground-truth signal.

% Sequence parameters.
flips = d2r(29.51); TR = 5.33e-3; Dur = 2.51e-3; Delta = 8058.48; n1B = 300; nMB = 300;
TBW = 2; B1rms = 4; dphi = 0;

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

%% Generate fitting dictionary.

%%% Load reconstruction dictionary here. %%%

n_increments = 100;
T1f = linspace(0.2,4,n_increments);
delta = 1;
M0s = linspace(0,0.2,n_increments);
T1s = 1;
T1d = linspace(0,8e-3,n_increments);
T2f = 84e-3;
T2s = 8.28e-6;
k = 55.2;
B0_var = 0;

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

for ii = 1:length(T2s)
    % Pre-compute lineshape for different T2s - can remove loop because T2s fixed.
    df_1B = [0 0 0]; df_2B = [0 0 Delta]; df_3B = [-Delta 0 Delta];
    [G_3B(ii,:),wloc_3B(ii)] = SuperLorentzian_LSint(T2s(ii),df_3B);
    [G_2B(ii,:),wloc_2B(ii)] = SuperLorentzian_LSint(T2s(ii),df_2B);
    [G_1B(ii,:),wloc_1B(ii)] = SuperLorentzian_LSint(T2s(ii),df_1B);
end

% Calculates signals for one sequence "period".
Mss_Total_noMT = zeros(size(Var_Params_noMT,1),5);
Mss_Total_MT = zeros(size(Var_Params_MT,1),5);
Mss_Total_ihMT = zeros(size(Var_Params_ihMT,1),5);

Dict_size = GetSize(Mss_Total_noMT) + GetSize(Mss_Total_MT) + GetSize(Mss_Total_ihMT);
disp(['Total dictionary size is: ',num2str(Dict_size/1e9),'GB.'])

% Signal generation.
delete(gcp('nocreate')); c = parcluster('local'); c.NumWorkers = 16; parpool(c, c.NumWorkers);

% Form M0 matrix for signal generation function.
M0_noMT_fit = [Var_Params_noMT(:,1) , Var_Params_noMT(:,2).*(1-Var_Params_noMT(:,3)) , Var_Params_noMT(:,2).*Var_Params_noMT(:,3)]; 
M0_MT_fit = [Var_Params_MT(:,1) , Var_Params_MT(:,2).*(1-Var_Params_MT(:,3)) , Var_Params_MT(:,2).*Var_Params_MT(:,3)]; 
M0_ihMT_fit = [Var_Params_ihMT(:,1) , Var_Params_ihMT(:,2).*(1-Var_Params_ihMT(:,3)) , Var_Params_ihMT(:,2).*Var_Params_ihMT(:,3)]; 

% No MT dictionary creation.
tic
parfor mm = 1:size(Var_Params_noMT,1)
    idx = find(T2s(:) == Var_Params_noMT(mm,9));
    tmp = Dictionary_function_CSS(flips,TR,Dur,Var_Params_noMT(mm,10),Var_Params_noMT(mm,5:7),Var_Params_noMT(mm,8:9),M0_noMT_fit(mm,:),Var_Params_noMT(mm,4),Delta,TBW,nMB,n1B,b1sqrd,G_1B(idx,:),G_2B(idx,:),G_3B(idx,:),wloc_1B(idx),wloc_2B(idx),wloc_3B(idx));
    Mss_Total_noMT(mm,:) = tmp*V_orig;
end
toc

% MT dictionary creation.
tic
parfor mm = 1:size(Var_Params_MT,1)
    idx = find(T2s(:) == Var_Params_MT(mm,9));
    tmp = Dictionary_function_CSS(flips,TR,Dur,Var_Params_MT(mm,10),Var_Params_MT(mm,5:7),Var_Params_MT(mm,8:9),M0_MT_fit(mm,:),Var_Params_MT(mm,4),Delta,TBW,nMB,n1B,b1sqrd,G_1B(idx,:),G_2B(idx,:),G_3B(idx,:),wloc_1B(idx),wloc_2B(idx),wloc_3B(idx));
    Mss_Total_MT(mm,:) = tmp*V_orig;
end
toc

% ihMT dictionary creation.
tic
parfor mm = 1:size(Var_Params_ihMT,1)
    idx = find(T2s(:) == Var_Params_ihMT(mm,9));
    tmp = Dictionary_function_CSS(flips,TR,Dur,Var_Params_ihMT(mm,10),Var_Params_ihMT(mm,5:7),Var_Params_ihMT(mm,8:9),M0_ihMT_fit(mm,:),Var_Params_ihMT(mm,4),Delta,TBW,nMB,n1B,b1sqrd,G_1B(idx,:),G_2B(idx,:),G_3B(idx,:),wloc_1B(idx),wloc_2B(idx),wloc_3B(idx));
    Mss_Total_ihMT(mm,:) = tmp*V_orig;
end
toc

Mss_AllSigs = [Mss_Total_noMT ; Mss_Total_MT ; Mss_Total_ihMT];
VP_AllSigs  = [Var_Params_noMT ; Var_Params_MT ; Var_Params_ihMT];

delete(gcp('nocreate'));
Mss_norm = conj(Mss_AllSigs./sum(abs(Mss_AllSigs).^2,2).^0.5); clear Mss_AllSigs
norm_atoms = Mss_norm; clear Mss_norm;

%% Test using default tissue parameters.

SNR = 40; nRealisations = 100;

M0b_WM = M0s(1:5:end);
T1D_WM = T1d(1:5:end)*1e3;

nsteps = length(M0b_WM);

param_fit = zeros(nsteps,nsteps,nRealisations,10); 
param_fit_norm = zeros(nsteps,nsteps,nRealisations,10); 
GT_vec = zeros(nsteps,nsteps,10);

for ii = 1:nsteps
    tic
    for jj = 1:nsteps
        
        % Ground-truth WM parameters.
        T1x_WM = [1090 1000 T1D_WM(jj)]*1e-3;
        T2x_WM = [84e-3, 8.28e-6];
        f_WM = 1;
        K_WM = 55.2;
        M0f_WM = 1-M0b_WM(ii);
        M0_WM = [M0f_WM M0b_WM(ii)*(1-f_WM) M0b_WM(ii)*f_WM];
        VP_WM = [M0_WM, T1x_WM, T2x_WM, K_WM];
        [GWM_3B,wlocWM_3B] = SuperLorentzian_LSint(T2x_WM(2),df_3B);
        [GWM_2B,wlocWM_2B] = SuperLorentzian_LSint(T2x_WM(2),df_2B);
        [GWM_1B,wlocWM_1B] = SuperLorentzian_LSint(T2x_WM(2),df_1B);
        
        x_WM = Dictionary_function_CSS(flips,TR,Dur,dphi,VP_WM(4:6),VP_WM(7:8),VP_WM(1:3),VP_WM(9),Delta,TBW,nMB,n1B,b1sqrd,GWM_1B,GWM_2B,GWM_3B,wlocWM_1B,wlocWM_2B,wlocWM_3B);
        GT_vec(ii,jj,:) = [M0f_WM, M0b_WM(ii), f_WM, K_WM, T1x_WM(1:2), T1D_WM(jj)*1e-3, T2x_WM, 0];
         
        for rr = 1:nRealisations
            
            %Add noise.
            xn_WM = zeros(1,1200); Sigma = mean(abs(x_WM))/SNR;
            for nn = 1:length(x_WM)
            xn_WM(nn) = x_WM(nn) + (normrnd(0,Sigma));
            end
            
            %Project to LR.
            xn_lr = xn_WM*V_orig;
            
            sig_norm = conj(xn_lr./sum(abs(xn_lr).^2,2).^0.5)';
            
            inner_product = norm_atoms * sig_norm;
            [~,max_idx] = max(inner_product,[],1);
            
            param_fit(ii,jj,rr,:) = VP_AllSigs(max_idx(:),:);
            param_fit_norm(ii,jj,rr,:) = ((squeeze(param_fit(ii,jj,rr,:))-squeeze(GT_vec(ii,jj,:)))./squeeze(GT_vec(ii,jj,:)))*100;
            
        end
        
    end
    toc
    disp(ii)
    
end

pf_mean = squeeze(mean(param_fit_norm,3));

figure(1); subplot(1,2,1)
imagesc(T1D_WM(4:end),M0b_WM(2:end),squeeze(pf_mean(2:end,4:end,2))); colormap(hot); 
tt = title('Relative Bias in Estimation of \it{f}'); tt.FontSize = 22;
set(gca,'FontSize',18); ylabel('\it{f}','FontSize',26); xlabel('T_{1D}^{s} (ms)','FontSize',20)
cb = colorbar; axis square; cb.Title.String = '(%)'; caxis([-5 5])
subplot(1,2,2)
imagesc(T1D_WM(4:end),M0b_WM(2:end),squeeze(pf_mean(2:end,4:end,7))); colormap(hot);
tt = title('Relative Bias in Estimation of T_{1D}^{s}'); tt.FontSize = 22;
set(gca,'FontSize',18); ylabel('\it{f}','FontSize',26); xlabel('T_{1D}^{s} (ms)','FontSize',20)
cb = colorbar; axis square; cb.Title.String = '(%)'; caxis([-5 5])
set(gcf,'Color','w','Name','Default Parameters');

figure(12); 
subplot(1,3,1)
imagesc(T1D_WM(4:end),M0b_WM(2:end),squeeze(mean(param_fit(2:end,4:end,:,5),3))); colormap(hot); 
tt = title('Estimates of T_{1}^{f}'); tt.FontSize = 22;
set(gca,'FontSize',18); ylabel('\it{f}','FontSize',26); xlabel('T_{1D}^{s} (ms)','FontSize',20)
cb = colorbar; axis square; cb.Title.String = '(s)'; caxis([min(T1f) max(T1f)])
subplot(1,3,2)
contourf(T1D_WM(4:end),M0b_WM(2:end),squeeze(mean(param_fit(2:end,4:end,:,2),3))); colormap(hot); 
tt = title('Estimates of \it{f}'); tt.FontSize = 22;
set(gca,'FontSize',18); ylabel('\it{f}','FontSize',26); xlabel('T_{1D}^{s} (ms)','FontSize',20)
cb = colorbar; axis square; caxis([min(M0b_WM) max(M0b_WM)])
subplot(1,3,3)
contourf(T1D_WM(4:end),M0b_WM(2:end),squeeze(mean(param_fit(2:end,4:end,:,7),3))*1e3); colormap(hot);
tt = title('Estimates of T_{1D}^{s}'); tt.FontSize = 22;
set(gca,'FontSize',18); ylabel('\it{f}','FontSize',26); xlabel('T_{1D}^{s} (ms)','FontSize',20)
cb = colorbar; axis square; cb.Title.String = '(ms)'; caxis([min(T1D_WM) max(T1D_WM)])
set(gcf,'Color','w','Name','Default Parameters');
