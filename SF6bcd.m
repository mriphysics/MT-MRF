%%% R2: Numerical simulation study. %%%

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

%% Test T1S

SNR = 40; nRealisations = 100;

M0b_WM = M0s(1:5:end);
T1D_WM = T1d(1:5:end)*1e3;

nsteps = length(M0b_WM); msteps = [600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400];

param_fit1 = zeros(length(msteps),nsteps,nsteps,nRealisations,10); 

for mm = 1:length(msteps)
    for ii = 1:nsteps
        tic
        for jj = 1:nsteps
            
            % Ground-truth WM parameters.
            T1x_WM = [1090 msteps(mm) T1D_WM(jj)]*1e-3;
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
                
                param_fit1(mm,ii,jj,rr,:) = VP_AllSigs(max_idx(:),:);
                
            end
            
        end
        toc
        disp(ii)
        
    end
end

%% Test T2F

SNR = 40; nRealisations = 100;

M0b_WM = M0s(1:5:end);
T1D_WM = T1d(1:5:end)*1e3;

nsteps = length(M0b_WM); msteps = [50e-3, 60e-3, 70e-3, 80e-3, 84e-3, 90e-3, 100e-3, 110e-3, 120e-3];

param_fit2 = zeros(length(msteps),nsteps,nsteps,nRealisations,10); 

for mm = 1:length(msteps)
    for ii = 1:nsteps
        tic
        for jj = 1:nsteps
            
            % Ground-truth WM parameters.
            T1x_WM = [1090 1000 T1D_WM(jj)]*1e-3;
            T2x_WM = [msteps(mm), 8.28e-6];
            f_WM = 1;
            K_WM = 55.2;
            M0f_WM = 1-M0b_WM(ii);
            M0_WM = [M0f_WM M0b_WM(ii)*(1-f_WM) M0b_WM(ii)*f_WM];
            VP_WM = [M0_WM, T1x_WM, T2x_WM, K_WM];
            [GWM_3B,wlocWM_3B] = SuperLorentzian_LSint(T2x_WM(2),df_3B);
            [GWM_2B,wlocWM_2B] = SuperLorentzian_LSint(T2x_WM(2),df_2B);
            [GWM_1B,wlocWM_1B] = SuperLorentzian_LSint(T2x_WM(2),df_1B);
            
            x_WM = Dictionary_function_CSS(flips,TR,Dur,dphi,VP_WM(4:6),VP_WM(7:8),VP_WM(1:3),VP_WM(9),Delta,TBW,nMB,n1B,b1sqrd,GWM_1B,GWM_2B,GWM_3B,wlocWM_1B,wlocWM_2B,wlocWM_3B);
            
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
                
                param_fit2(mm,ii,jj,rr,:) = VP_AllSigs(max_idx(:),:);
                
            end
            
        end
        toc
        disp(ii)
        
    end
end

%% Test T2S

SNR = 40; nRealisations = 100;

M0b_WM = M0s(1:5:end);
T1D_WM = T1d(1:5:end)*1e3;

nsteps = length(M0b_WM); msteps = [6e-6, 7e-6, 8e-6, 8.28e-6, 9e-6, 10e-6, 11e-6, 12e-6, 13e-6];

param_fit3 = zeros(length(msteps),nsteps,nsteps,nRealisations,10); 

for mm = 1:length(msteps)
    for ii = 1:nsteps
        tic
        for jj = 1:nsteps
            
            % Ground-truth WM parameters.
            T1x_WM = [1090 1000 T1D_WM(jj)]*1e-3;
            T2x_WM = [84e-3, msteps(mm)];
            f_WM = 1;
            K_WM = 55.2;
            M0f_WM = 1-M0b_WM(ii);
            M0_WM = [M0f_WM M0b_WM(ii)*(1-f_WM) M0b_WM(ii)*f_WM];
            VP_WM = [M0_WM, T1x_WM, T2x_WM, K_WM];
            [GWM_3B,wlocWM_3B] = SuperLorentzian_LSint(T2x_WM(2),df_3B);
            [GWM_2B,wlocWM_2B] = SuperLorentzian_LSint(T2x_WM(2),df_2B);
            [GWM_1B,wlocWM_1B] = SuperLorentzian_LSint(T2x_WM(2),df_1B);
            
            x_WM = Dictionary_function_CSS(flips,TR,Dur,dphi,VP_WM(4:6),VP_WM(7:8),VP_WM(1:3),VP_WM(9),Delta,TBW,nMB,n1B,b1sqrd,GWM_1B,GWM_2B,GWM_3B,wlocWM_1B,wlocWM_2B,wlocWM_3B);
            
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
                
                param_fit3(mm,ii,jj,rr,:) = VP_AllSigs(max_idx(:),:);
                
            end
            
        end
        toc
        disp(ii)
        
    end
end

%% Test K

SNR = 40; nRealisations = 100;

M0b_WM = M0s(1:5:end);
T1D_WM = T1d(1:5:end)*1e3;

nsteps = length(M0b_WM); msteps = [45, 50, 55, 55.2, 60, 65, 70, 75, 80];

param_fit4 = zeros(length(msteps),nsteps,nsteps,nRealisations,10); 

for mm = 1:length(msteps)
    for ii = 1:nsteps
        tic
        for jj = 1:nsteps
            
            % Ground-truth WM parameters.
            T1x_WM = [1090 1000 T1D_WM(jj)]*1e-3;
            T2x_WM = [84e-3, 8.28e-6];
            f_WM = 1;
            K_WM = msteps(mm);
            M0f_WM = 1-M0b_WM(ii);
            M0_WM = [M0f_WM M0b_WM(ii)*(1-f_WM) M0b_WM(ii)*f_WM];
            VP_WM = [M0_WM, T1x_WM, T2x_WM, K_WM];
            [GWM_3B,wlocWM_3B] = SuperLorentzian_LSint(T2x_WM(2),df_3B);
            [GWM_2B,wlocWM_2B] = SuperLorentzian_LSint(T2x_WM(2),df_2B);
            [GWM_1B,wlocWM_1B] = SuperLorentzian_LSint(T2x_WM(2),df_1B);
            
            x_WM = Dictionary_function_CSS(flips,TR,Dur,dphi,VP_WM(4:6),VP_WM(7:8),VP_WM(1:3),VP_WM(9),Delta,TBW,nMB,n1B,b1sqrd,GWM_1B,GWM_2B,GWM_3B,wlocWM_1B,wlocWM_2B,wlocWM_3B);
            
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
                
                param_fit4(mm,ii,jj,rr,:) = VP_AllSigs(max_idx(:),:);
                
            end
            
        end
        toc
        disp(ii)
        
    end
end

%% Generate heat maps.

str = {'T_{1Z}^{s} = 600ms','T_{1Z}^{s} = 700ms',...
       'T_{1Z}^{s} = 800ms','T_{1Z}^{s} = 900ms',...
       'T_{1Z}^{s} = 1000ms (GT)','T_{1Z}^{s} = 1100ms',...
       'T_{1Z}^{s} = 1200ms','T_{1Z}^{s} = 1300ms','T_{1Z}^{s} = 1400ms'};

for mm = 1:length(msteps)
    
    figure('units','normalized','outerposition',[0 0 1 1])
    subplot(1,3,1)
    imagesc(T1D_WM(4:end),M0b_WM(2:end),squeeze(mean(param_fit1(mm,2:end,4:end,:,5),4))); colormap(hot);
    tt = title('T_{1}^{f} estimates'); tt.FontSize = 22;
    set(gca,'FontSize',18); ylabel('\it{f}','FontSize',26); xlabel('T_{1D}^{s} (ms)','FontSize',20)
    cb = colorbar; cb.Title.String = '(s)'; axis square; caxis([min(T1f) max(T1f)])
    subplot(1,3,2)
    contourf(T1D_WM(4:end),M0b_WM(2:end),squeeze(mean(param_fit1(mm,2:end,4:end,:,2),4))); colormap(hot);
    tt = title('{\itf} estimates'); tt.FontSize = 22;
    set(gca,'FontSize',18); ylabel('\it{f}','FontSize',26); xlabel('T_{1D}^{s} (ms)','FontSize',20)
    cb = colorbar; axis square; caxis([min(M0b_WM) max(M0b_WM)])
    subplot(1,3,3)
    contourf(T1D_WM(4:end),M0b_WM(2:end),squeeze(mean(param_fit1(mm,2:end,4:end,:,7),4))*1e3); colormap(hot);
    tt = title('T_{1D}^{s} estimates'); tt.FontSize = 22;
    set(gca,'FontSize',18); ylabel('\it{f}','FontSize',26); xlabel('T_{1D}^{s} (ms)','FontSize',20)
    cb = colorbar; axis square; cb.Title.String = '(ms)'; caxis([min(T1D_WM) max(T1D_WM)])
    set(gcf,'Color','w');
    annotation('textbox','String',str(mm),'FontSize',28,'Position',...
        [0.01 0.922804202992474 0.243716599866399 0.0558482599773443],'LineStyle','none','FontWeight','bold',...
        'FitBoxToText','off','HorizontalAlignment','left')
    set(gcf,'Position',[0 0.283 0.998 0.617])
    saveas(gcf,sprintf('varT1S%d.png',mm))

end

str = {'T_{2}^{f} = 50ms','T_{2}^{f} = 60ms',...
       'T_{2}^{f} = 70ms','T_{2}^{f} = 80ms',...
       'T_{2}^{f} = 84ms (GT)','T_{2}^{f} = 90ms',...
       'T_{2}^{f} = 100ms','T_{2}^{f} = 110ms','T_{2}^{f} = 120ms'};
   
for mm = 1:length(msteps)
    
    figure('units','normalized','outerposition',[0 0 1 1])
    subplot(1,3,1)
    imagesc(T1D_WM(4:end),M0b_WM(2:end),squeeze(mean(param_fit2(mm,2:end,4:end,:,5),4))); colormap(hot);
    tt = title('T_{1}^{f} estimates'); tt.FontSize = 22;
    set(gca,'FontSize',18); ylabel('\it{f}','FontSize',26); xlabel('T_{1D}^{s} (ms)','FontSize',20)
    cb = colorbar; cb.Title.String = '(s)'; axis square; caxis([min(T1f) max(T1f)])
    subplot(1,3,2)
    contourf(T1D_WM(4:end),M0b_WM(2:end),squeeze(mean(param_fit2(mm,2:end,4:end,:,2),4))); colormap(hot);
    tt = title('{\itf} estimates'); tt.FontSize = 22;
    set(gca,'FontSize',18); ylabel('\it{f}','FontSize',26); xlabel('T_{1D}^{s} (ms)','FontSize',20)
    cb = colorbar; axis square; caxis([min(M0b_WM) max(M0b_WM)])
    subplot(1,3,3)
    contourf(T1D_WM(4:end),M0b_WM(2:end),squeeze(mean(param_fit2(mm,2:end,4:end,:,7),4))*1e3); colormap(hot);
    tt = title('T_{1D}^{s} estimates'); tt.FontSize = 22;
    set(gca,'FontSize',18); ylabel('\it{f}','FontSize',26); xlabel('T_{1D}^{s} (ms)','FontSize',20)
    cb = colorbar; axis square; cb.Title.String = '(ms)'; caxis([min(T1D_WM) max(T1D_WM)])
    set(gcf,'Color','w');
    annotation('textbox','String',str(mm),'FontSize',28,'Position',...
        [0.01 0.922804202992474 0.243716599866399 0.0558482599773443],'LineStyle','none','FontWeight','bold',...
        'FitBoxToText','off','HorizontalAlignment','left')
    set(gcf,'Position',[0 0.283 0.998 0.617])   
    saveas(gcf,sprintf('varT2F%d.png',mm))

end

str = {'T_{2}^{s} = 6\mus','T_{2}^{s} = 7\mus',...
       'T_{2}^{s} = 8\mus','T_{2}^{s} = 8.28\mus (GT)',...
       'T_{2}^{s} = 9\mus','T_{2}^{s} = 10\mus',...
       'T_{2}^{s} = 11\mus','T_{2}^{s} = 12\mus','T_{2}^{s} = 13\mus'};

for mm = 1:length(msteps)
    
    figure('units','normalized','outerposition',[0 0 1 1])
    subplot(1,3,1)
    imagesc(T1D_WM(4:end),M0b_WM(2:end),squeeze(mean(param_fit3(mm,2:end,4:end,:,5),4))); colormap(hot);
    tt = title('T_{1}^{f} estimates'); tt.FontSize = 22;
    set(gca,'FontSize',18); ylabel('\it{f}','FontSize',26); xlabel('T_{1D}^{s} (ms)','FontSize',20)
    cb = colorbar; cb.Title.String = '(s)'; axis square; caxis([min(T1f) max(T1f)])
    subplot(1,3,2)
    contourf(T1D_WM(4:end),M0b_WM(2:end),squeeze(mean(param_fit3(mm,2:end,4:end,:,2),4))); colormap(hot);
    tt = title('{\itf} estimates'); tt.FontSize = 22;
    set(gca,'FontSize',18); ylabel('\it{f}','FontSize',26); xlabel('T_{1D}^{s} (ms)','FontSize',20)
    cb = colorbar; axis square; caxis([min(M0b_WM) max(M0b_WM)])
    subplot(1,3,3)
    contourf(T1D_WM(4:end),M0b_WM(2:end),squeeze(mean(param_fit3(mm,2:end,4:end,:,7),4))*1e3); colormap(hot);
    tt = title('T_{1D}^{s} estimates'); tt.FontSize = 22;
    set(gca,'FontSize',18); ylabel('\it{f}','FontSize',26); xlabel('T_{1D}^{s} (ms)','FontSize',20)
    cb = colorbar; axis square; cb.Title.String = '(ms)'; caxis([min(T1D_WM) max(T1D_WM)])
    set(gcf,'Color','w');
    annotation('textbox','String',str(mm),'FontSize',28,'Position',...
        [0.01 0.922804202992474 0.243716599866399 0.0558482599773443],'LineStyle','none','FontWeight','bold',...
        'FitBoxToText','off','HorizontalAlignment','left')
    set(gcf,'Position',[0 0.283 0.998 0.617])   
    saveas(gcf,sprintf('varT2S%d.png',mm))

end

str = {'{\itk} = 45s^{-1}','{\itk} = 50s^{-1}',...
       '{\itk} = 55s^{-1}','{\itk} = 55.2s^{-1} (GT)',...
       '{\itk} = 60s^{-1}','{\itk} = 65s^{-1}',...
       '{\itk} = 70s^{-1}','{\itk} = 75s^{-1}','{\itk} = 80s^{-1}'};

for mm = 1:length(msteps)
    
    figure('units','normalized','outerposition',[0 0 1 1])
    subplot(1,3,1)
    imagesc(T1D_WM(4:end),M0b_WM(2:end),squeeze(mean(param_fit4(mm,2:end,4:end,:,5),4))); colormap(hot);
    tt = title('T_{1}^{f} estimates'); tt.FontSize = 22;
    set(gca,'FontSize',18); ylabel('\it{f}','FontSize',26); xlabel('T_{1D}^{s} (ms)','FontSize',20)
    cb = colorbar; cb.Title.String = '(s)'; axis square; caxis([min(T1f) max(T1f)])
    subplot(1,3,2)
    contourf(T1D_WM(4:end),M0b_WM(2:end),squeeze(mean(param_fit4(mm,2:end,4:end,:,2),4))); colormap(hot);
    tt = title('{\itf} estimates'); tt.FontSize = 22;
    set(gca,'FontSize',18); ylabel('\it{f}','FontSize',26); xlabel('T_{1D}^{s} (ms)','FontSize',20)
    cb = colorbar; axis square; caxis([min(M0b_WM) max(M0b_WM)])
    subplot(1,3,3)
    contourf(T1D_WM(4:end),M0b_WM(2:end),squeeze(mean(param_fit4(mm,2:end,4:end,:,7),4))*1e3); colormap(hot);
    tt = title('T_{1D}^{s} estimates'); tt.FontSize = 22;
    set(gca,'FontSize',18); ylabel('\it{f}','FontSize',26); xlabel('T_{1D}^{s} (ms)','FontSize',20)
    cb = colorbar; axis square; cb.Title.String = '(ms)'; caxis([min(T1D_WM) max(T1D_WM)])
    set(gcf,'Color','w');
    annotation('textbox','String',str(mm),'FontSize',28,'Position',...
        [0.01 0.922804202992474 0.243716599866399 0.0558482599773443],'LineStyle','none','FontWeight','bold',...
        'FitBoxToText','off','HorizontalAlignment','left')
    set(gcf,'Position',[0 0.283 0.998 0.617])   
    saveas(gcf,sprintf('vark%d.png',mm))

end

close all;

%% Generate scatter plots.

M0b_WM = M0s(1:5:end);
T1D_WM = T1d(1:5:end)*1e3;
T1F_GT = 1.090;

ss_grid = repmat(M0b_WM.',[1 length(M0b_WM)]); ss_grid = ss_grid(2:end,4:end);
T1d_grid = repmat((T1D_WM).',[1 length(M0b_WM)]).'; T1d_grid = T1d_grid(2:end,4:end);

ss_grid_all = zeros(size(param_fit1,1),size(param_fit1,2)-1,size(param_fit1,3)-3,size(param_fit1,4)); 
T1d_grid_all = zeros(size(param_fit1,1),size(param_fit1,2)-1,size(param_fit1,3)-3,size(param_fit1,4)); 

for ii = 1:size(param_fit1,1)
    for jj = 1:size(param_fit1,4)
        ss_grid_all(ii,:,:,jj) = ss_grid;
        T1d_grid_all(ii,:,:,jj) = T1d_grid;
    end
end

T1f_grid_all = T1F_GT*ones(size(param_fit1,1),size(param_fit1,2)-1,size(param_fit1,3)-3,size(param_fit1,4)); 

param_fit1_f = reshape(param_fit1(:,2:end,4:end,:,2)-ss_grid_all,[size(ss_grid_all,1) size(ss_grid_all,2)*size(ss_grid_all,3)*size(ss_grid_all,4)]);
param_fit1_T1d = reshape(param_fit1(:,2:end,4:end,:,7)*1e3-T1d_grid_all,[size(ss_grid_all,1) size(ss_grid_all,2)*size(ss_grid_all,3)*size(ss_grid_all,4)]);
param_fit1_T1f = reshape(param_fit1(:,2:end,4:end,:,5)-T1f_grid_all,[size(ss_grid_all,1) size(ss_grid_all,2)*size(ss_grid_all,3)*size(ss_grid_all,4)]);

param_fit2_f = reshape(param_fit2(:,2:end,4:end,:,2)-ss_grid_all,[size(ss_grid_all,1) size(ss_grid_all,2)*size(ss_grid_all,3)*size(ss_grid_all,4)]);
param_fit2_T1d = reshape(param_fit2(:,2:end,4:end,:,7)*1e3-T1d_grid_all,[size(ss_grid_all,1) size(ss_grid_all,2)*size(ss_grid_all,3)*size(ss_grid_all,4)]);
param_fit2_T1f = reshape(param_fit2(:,2:end,4:end,:,5)-T1f_grid_all,[size(ss_grid_all,1) size(ss_grid_all,2)*size(ss_grid_all,3)*size(ss_grid_all,4)]);

param_fit3_f = reshape(param_fit3(:,2:end,4:end,:,2)-ss_grid_all,[size(ss_grid_all,1) size(ss_grid_all,2)*size(ss_grid_all,3)*size(ss_grid_all,4)]);
param_fit3_T1d = reshape(param_fit3(:,2:end,4:end,:,7)*1e3-T1d_grid_all,[size(ss_grid_all,1) size(ss_grid_all,2)*size(ss_grid_all,3)*size(ss_grid_all,4)]);
param_fit3_T1f = reshape(param_fit3(:,2:end,4:end,:,5)-T1f_grid_all,[size(ss_grid_all,1) size(ss_grid_all,2)*size(ss_grid_all,3)*size(ss_grid_all,4)]);

param_fit4_f = reshape(param_fit4(:,2:end,4:end,:,2)-ss_grid_all,[size(ss_grid_all,1) size(ss_grid_all,2)*size(ss_grid_all,3)*size(ss_grid_all,4)]);
param_fit4_T1d = reshape(param_fit4(:,2:end,4:end,:,7)*1e3-T1d_grid_all,[size(ss_grid_all,1) size(ss_grid_all,2)*size(ss_grid_all,3)*size(ss_grid_all,4)]);
param_fit4_T1f = reshape(param_fit4(:,2:end,4:end,:,5)-T1f_grid_all,[size(ss_grid_all,1) size(ss_grid_all,2)*size(ss_grid_all,3)*size(ss_grid_all,4)]);

cm = colormap(linspecer(3));
figure(1)
subplot(3,1,1); 
violinplot(param_fit1_T1f.',cellstr({'600','700','800','900','1000','1100','1200','1300','1400'}),'BoxColor',cm(1,:),'MedianColor',cm(1,:),'EdgeColor',cm(1,:),'ViolinColor',cm(1,:),'ViolinAlpha',0.6);
set(gca,'FontSize',18); ylabel('T_{1}^{f} error (s)','FontSize',20); xlabel('T_{1Z}^{s} (ms)','FontSize',20); grid on; grid minor;
ylim([-0.4 0.4])
subplot(3,1,2); 
violinplot(param_fit1_f.',cellstr({'600','700','800','900','1000','1100','1200','1300','1400'}),'BoxColor',cm(2,:),'MedianColor',cm(2,:),'EdgeColor',cm(2,:),'ViolinColor',cm(2,:),'ViolinAlpha',0.6);
set(gca,'FontSize',18); ylabel('{\itf} error','FontSize',20); xlabel('T_{1Z}^{s} (ms)','FontSize',20); grid on; grid minor;
ylim([-0.06 0.03])
subplot(3,1,3); 
violinplot(param_fit1_T1d.',cellstr({'600','700','800','900','1000','1100','1200','1300','1400'}),'BoxColor',cm(3,:),'MedianColor',cm(3,:),'EdgeColor',cm(3,:),'ViolinColor',cm(3,:),'ViolinAlpha',0.6);
set(gca,'FontSize',18); ylabel('T_{1D}^{s} error (ms)','FontSize',20); xlabel('T_{1Z}^{s} (ms)','FontSize',20); grid on; grid minor;
ylim([-5 5])
set(gcf,'units','normalized','outerposition',[0 0 1 1]); saveas(gcf,sprintf('violinT1s.png'))

figure(2)
subplot(3,1,1); 
violinplot(param_fit2_T1f.',cellstr({'50','60','70','80','84','90','100','110','120'}),'BoxColor',cm(1,:),'MedianColor',cm(1,:),'EdgeColor',cm(1,:),'ViolinColor',cm(1,:),'ViolinAlpha',0.6);
set(gca,'FontSize',18); ylabel('T_{1}^{f} error (s)','FontSize',20); xlabel('T_{2}^{f} (ms)','FontSize',20); grid on; grid minor;
ylim([-0.4 0.4])
subplot(3,1,2); 
violinplot(param_fit2_f.',cellstr({'50','60','70','80','84','90','100','110','120'}),'BoxColor',cm(2,:),'MedianColor',cm(2,:),'EdgeColor',cm(2,:),'ViolinColor',cm(2,:),'ViolinAlpha',0.6);
set(gca,'FontSize',18); ylabel('{\itf} error','FontSize',20); xlabel('T_{2}^{f} (ms)','FontSize',20); grid on; grid minor;
ylim([-0.06 0.03])
subplot(3,1,3); 
violinplot(param_fit2_T1d.',cellstr({'50','60','70','80','84','90','100','110','120'}),'BoxColor',cm(3,:),'MedianColor',cm(3,:),'EdgeColor',cm(3,:),'ViolinColor',cm(3,:),'ViolinAlpha',0.6);
set(gca,'FontSize',18); ylabel('T_{1D}^{s} error (ms)','FontSize',20); xlabel('T_{2}^{f} (ms)','FontSize',20); grid on; grid minor;
ylim([-5 5])
set(gcf,'units','normalized','outerposition',[0 0 1 1]); saveas(gcf,sprintf('violinT2f.png'))

figure(3)
subplot(3,1,1); 
violinplot(param_fit3_T1f.',cellstr({'6','7','8','8.28','9','10','11','12','13'}),'BoxColor',cm(1,:),'MedianColor',cm(1,:),'EdgeColor',cm(1,:),'ViolinColor',cm(1,:),'ViolinAlpha',0.6);
set(gca,'FontSize',18); ylabel('T_{1}^{f} error (s)','FontSize',20); xlabel('T_{2}^{s} (\mus)','FontSize',20); grid on; grid minor;
ylim([-0.4 0.4])
subplot(3,1,2); 
violinplot(param_fit3_f.',cellstr({'6','7','8','8.28','9','10','11','12','13'}),'BoxColor',cm(2,:),'MedianColor',cm(2,:),'EdgeColor',cm(2,:),'ViolinColor',cm(2,:),'ViolinAlpha',0.6);
set(gca,'FontSize',18); ylabel('{\itf} error','FontSize',20); xlabel('T_{2}^{s} (\mus)','FontSize',20); grid on; grid minor;
ylim([-0.06 0.03])
subplot(3,1,3); 
violinplot(param_fit3_T1d.',cellstr({'6','7','8','8.28','9','10','11','12','13'}),'BoxColor',cm(3,:),'MedianColor',cm(3,:),'EdgeColor',cm(3,:),'ViolinColor',cm(3,:),'ViolinAlpha',0.6);
set(gca,'FontSize',18); ylabel('T_{1D}^{s} error (ms)','FontSize',20); xlabel('T_{2}^{s} (\mus)','FontSize',20); grid on; grid minor;
ylim([-5 5])
set(gcf,'units','normalized','outerposition',[0 0 1 1]); saveas(gcf,sprintf('violinT2s.png'))

figure(4)
subplot(3,1,1); 
violinplot(param_fit4_T1f.',cellstr({'45','50','55','55.2','60','65','70','75','80'}),'BoxColor',cm(1,:),'MedianColor',cm(1,:),'EdgeColor',cm(1,:),'ViolinColor',cm(1,:),'ViolinAlpha',0.6);
set(gca,'FontSize',18); ylabel('T_{1}^{f} error (s)','FontSize',20); xlabel('{\itk} (s^{-1})','FontSize',20); grid on; grid minor;
ylim([-0.4 0.4])
subplot(3,1,2); 
violinplot(param_fit4_f.',cellstr({'45','50','55','55.2','60','65','70','75','80'}),'BoxColor',cm(2,:),'MedianColor',cm(2,:),'EdgeColor',cm(2,:),'ViolinColor',cm(2,:),'ViolinAlpha',0.6);
set(gca,'FontSize',18); ylabel('{\itf} error','FontSize',20); xlabel('{\itk} (s^{-1})','FontSize',20); grid on; grid minor;
ylim([-0.06 0.03])
subplot(3,1,3); 
violinplot(param_fit4_T1d.',cellstr({'45','50','55','55.2','60','65','70','75','80'}),'BoxColor',cm(3,:),'MedianColor',cm(3,:),'EdgeColor',cm(3,:),'ViolinColor',cm(3,:),'ViolinAlpha',0.6);
set(gca,'FontSize',18); ylabel('T_{1D}^{s} error (ms)','FontSize',20); xlabel('{\itk} (s^{-1})','FontSize',20); grid on; grid minor;
ylim([-5 5])
set(gcf,'units','normalized','outerposition',[0 0 1 1]); saveas(gcf,sprintf('violink.png'))
