%% Example script for dictionary matching in SVD space. Daniel West 2020

% Load data.
load('/home/dw16/TGR_files/hdp3d_ex.mat'); % In-vivo data.
%load('/home/dw16/TGR_files/hdp3d_ext.mat'); % Phantom data.

% Load dictionary.
load('/home/dw16/TGR_files/Phantom_Dict_ex2.mat'); V_orig = V(:,1:5);

if exist('hdp3d_ex','var') == 1
    hdp3d = hdp3d_ex; clear hdp3d_ex;
end

Mss_norm = conj(Mss_AllSigs./sum(abs(Mss_AllSigs).^2,2).^0.5); clear Mss_AllSigs
norm_atoms = Mss_norm*V_orig; clear Mss_norm;

% Project to time-domain.
MRF_sing = hdp3d(84:192,95:240,7:116,:); % In-vivo data.
%MRF_sing = squeeze(hdp3d(63:90,50:105,23:74,:)); % Phantom data.
nx = size(MRF_sing,1); ny = size(MRF_sing,2); nz = size(MRF_sing,3); ns = size(MRF_sing,4); nt = size(V_orig,1);
mrf_aux = reshape(MRF_sing,[nx*ny*nz ns]); mrf_aux = mrf_aux.';
MRF_LRI5 = (V_orig*mrf_aux).'; MRF_LRI5 = reshape(MRF_LRI5,[nx ny nz nt]);

% Take out global phase.
phase = repmat(angle(sum(MRF_LRI5,4)),[1 1 1 nt]);
MRF_LRI5 = MRF_LRI5.*exp(-1i*phase);

% Go back to SVD domain.
mrf_aux = reshape(MRF_LRI5,[nx*ny*nz nt]);
MRF_LRI5 = (mrf_aux*V_orig); MRF_LRI5 = reshape(MRF_LRI5,[nx ny nz ns]);

data_dim = size(MRF_LRI5);
MRF_LRI5_reshaped = permute(MRF_LRI5,[4 1 2 3]);
MRF_LRI5_reshaped = reshape(MRF_LRI5_reshaped,[data_dim(4), data_dim(1)*data_dim(2)*data_dim(3)]);

npos = data_dim(1)*data_dim(2)*data_dim(3); max_idx = zeros(npos,1);

sig_norm = conj(MRF_LRI5_reshaped./sum(abs(MRF_LRI5_reshaped).^2,1).^0.5);

nchunks = ceil(size(norm_atoms,1)*size(sig_norm,2)/(20*1024*1024*128));

if nchunks == 1
    tic
    inner_product = norm_atoms * sig_norm;
    toc
    % Find maximum of each column.
    [~,max_idx] = max(inner_product,[],1);
else
    size_chunk = floor(npos/nchunks);
    for cc=1:nchunks
        if cc == nchunks
            aux = (cc-1)*size_chunk+1:npos;
        else
            aux = (cc-1)*size_chunk+1:cc*size_chunk;
        end
        inner_product = norm_atoms * sig_norm(:,aux);
        
        % Find maximum of each column.
        [~,max_idx(aux)] = max(inner_product,[],1);
        
        clear inner_product
        
        pause(5)
        
    end   
end

% Find dictionary entry.
param_fit = VP_AllSigs(max_idx(:),:);
params_reshaped = reshape(param_fit,[data_dim(1), data_dim(2), data_dim(3) 10]);        

%% Create final parameter maps.

clearvars -except norm_atoms VP_AllSigs params_reshaped

load('/home/dw16/TGR_files/hdp3d_ex.mat'); load('/home/dw16/TGR_files/Phantom_Dict_ex2.mat'); V_orig = V(:,1:5); % In-vivo data.
%load('/home/dw16/TGR_files/hdp3d_ext.mat'); load('/home/dw16/TGR_files/Phantom_Dict_ex2.mat'); V_orig = V(:,1:5); % Phantom data.

hdp3d = hdp3d_ex(84:192,95:240,7:116,:); % In-vivo data.
%hdp3d = hdp3d(63:90,50:105,23:74,:); % Phantom data.
mrf_aux = reshape(hdp3d,[size(hdp3d,1)*size(hdp3d,2)*size(hdp3d,3) size(hdp3d,4)]); mrf_aux = mrf_aux';
MRF_LRI5 = (V_orig*mrf_aux)'; MRF_LRI5 = reshape(MRF_LRI5,[size(hdp3d,1) size(hdp3d,2) size(hdp3d,3) size(MRF_LRI5,2)]);
mask = squeeze(abs(MRF_LRI5(:,:,:,600))) > 1.5;

M0s_map = squeeze(params_reshaped(:,:,:,2)).*mask; delta_map = squeeze(params_reshaped(:,:,:,7)).*mask; T1_map = squeeze(params_reshaped(:,:,:,5)).*mask;
%K_map = squeeze(params_reshaped(:,:,:,4)).*mask; % Initialise if phantom data is being analysed.

%% Load and visualise segmentation for in-vivo analysis.

% Convert M0b_map to NIFTI and use FSL BET + FAST. Loading output:
niftiwrite(M0s_map,'M0s_map.nii');

CSF_mask = load_nii('/home/dw16/TGR_files/M0s_map_BET_pve_0.nii.gz'); CSF_mask = CSF_mask.img;
GM_mask = load_nii('/home/dw16/TGR_files/M0s_map_BET_pve_1.nii.gz'); GM_mask = GM_mask.img;
WM_mask = load_nii('/home/dw16/TGR_files/M0s_map_BET_pve_2.nii.gz'); WM_mask = WM_mask.img;

figure(1); 
subplot(3,3,1); imagesc(squeeze(CSF_mask(:,:,36))); axis image; set(gca,'YTickLabels',[],'visible','off'); set(gca,'XTickLabels',[],'visible','off');
subplot(3,3,4); imagesc(squeeze(GM_mask(:,:,36))); axis image; set(gca,'YTickLabels',[],'visible','off'); set(gca,'XTickLabels',[],'visible','off');
subplot(3,3,7); imagesc(squeeze(WM_mask(:,:,36))); axis image; set(gca,'YTickLabels',[],'visible','off'); set(gca,'XTickLabels',[],'visible','off');
subplot(3,3,2); imagesc(squeeze(CSF_mask(67,:,:))); axis image; set(gca,'YTickLabels',[],'visible','off'); set(gca,'XTickLabels',[],'visible','off');
subplot(3,3,5); imagesc(squeeze(GM_mask(67,:,:))); axis image; set(gca,'YTickLabels',[],'visible','off'); set(gca,'XTickLabels',[],'visible','off');
subplot(3,3,8); imagesc(squeeze(WM_mask(67,:,:))); axis image; set(gca,'YTickLabels',[],'visible','off'); set(gca,'XTickLabels',[],'visible','off');
subplot(3,3,3); imagesc(squeeze(CSF_mask(:,71,:))); axis image; set(gca,'YTickLabels',[],'visible','off'); set(gca,'XTickLabels',[],'visible','off');
subplot(3,3,6); imagesc(squeeze(GM_mask(:,71,:))); axis image; set(gca,'YTickLabels',[],'visible','off'); set(gca,'XTickLabels',[],'visible','off');
subplot(3,3,9); imagesc(squeeze(WM_mask(:,71,:))); axis image; set(gca,'YTickLabels',[],'visible','off'); set(gca,'XTickLabels',[],'visible','off');
set(gcf,'color','k','InvertHardCopy','off'); colormap(magma());

%% Create joint histograms.

Slicex = 20:90; Slicey = 20:125; Slicez = 15:95;
WM_idx = find(WM_mask(Slicex,Slicey,Slicez) == 1); [WMx, WMy, WMz] = ind2sub(size(WM_mask(Slicex,Slicey,Slicez)),WM_idx);
GM_idx = find(GM_mask(Slicex,Slicey,Slicez) == 1); [GMx, GMy, GMz] = ind2sub(size(GM_mask(Slicex,Slicey,Slicez)),GM_idx);

delta_WM = zeros(length(WMx),1); M0s_WM = zeros(length(WMx),1);
for ii = 1:length(WMx)
    delta_WM(ii) = squeeze(delta_map(WMx(ii)+min(Slicex)-1,min(Slicey)+WMy(ii)-1,min(Slicez)+WMz(ii)-1));
    M0s_WM(ii) = squeeze(M0s_map(WMx(ii)+min(Slicex)-1,min(Slicey)+WMy(ii)-1,min(Slicez)+WMz(ii)-1));
end

delta_GM = zeros(length(GMx),1); M0s_GM = zeros(length(GMx),1);
for ii = 1:length(GMx)
    delta_GM(ii) = squeeze(delta_map(GMx(ii)+min(Slicex)-1,min(Slicey)+GMy(ii)-1,min(Slicez)+GMz(ii)-1));
    M0s_GM(ii) = squeeze(M0s_map(GMx(ii)+min(Slicex)-1,min(Slicey)+GMy(ii)-1,min(Slicez)+GMz(ii)-1));
end

nbins = 22; delta_vector = linspace(1e-3,7e-3,nbins+5); M0s_vector = linspace(0.03,0.11,nbins);
hist_example = histcounts2([delta_WM(:);delta_GM(:)], [M0s_WM(:);M0s_GM(:)], delta_vector, M0s_vector);

fig = figure;
subplot('Position', [0.3, 0.35, 0.6, 0.6]);
im = imagesc(delta_vector*1e3,M0s_vector,hist_example.'); colormap(plasma); cb = colorbar; caxis([0 2500]);
set(gca,'FontSize',20); xlabel('T_{1D}^{s} (ms)','FontSize',20); ylabel('\it{f}','FontSize',26); cb.FontSize = 20; 
title(cb,'Count','FontSize',20,'Position',[51.49999999999999,412.7500131955103,0]);
cb.Ticks = 0:500:2500; cb.TickLabels = {0:500:2500}; set(gcf,'Units','Normalized');
im.Parent.YDir = 'normal'; pos = get(gcf,'Position');

subplot('Position', [0.43, 0.05, 0.45, 0.15]);
histogram(delta_WM, delta_vector,'FaceColor',[0.6 0.6 0.6],'EdgeColor',[0.6 0.6 0.6]); hold on
histogram(delta_GM, delta_vector,'FaceColor',[0.1 0.1 0.1],'EdgeColor',[0.1 0.1 0.1]);
camroll(180); set(gca,'FontSize',20); grid on; grid minor;
set(gca,'XDir','reverse','Position',[0.273 0.060 0.625 0.185]);
set(gca,'XTickLabels',[],'Visible','off'); ylim([0 1.2e4])

subplot('Position', [0, 0.4, 0.15, 0.55]);
histogram(M0s_WM, M0s_vector,'FaceColor',[0.6 0.6 0.6],'EdgeColor',[0.6 0.6 0.6]); hold on
histogram(M0s_GM, M0s_vector,'FaceColor',[0.1 0.1 0.1],'EdgeColor',[0.1 0.1 0.1]);
ll = legend('WM','GM'); ll.FontSize = 20; ll.Position = [0.07996527855388,0.875504041349496,0.052604165890565,0.072181241404269]; legend boxoff;
camroll(90); set(gca,'FontSize',20); grid on; grid minor; ylim([0 3e4])
set(gca,'Position',[0.015 0.321 0.220 0.659]);
set(gca,'XTickLabels',[],'Visible','off'); 

%% Plot final maps. Commented out sections are for phantom data.

Sigma = 0.7;  %Slice = 28;
figure(2)
subplot(3,3,2); imagesc(imgaussfilt(squeeze(T1_map(67,:,:)),Sigma)); axis image; 
colormap('magma'); cb = colorbar; cb.FontSize = 20; cb.Color = 'w'; title(cb,'s','FontSize',20,'Color','w'); caxis([0 2.5]);
set(gca,'YTickLabels',[],'visible','off'); set(gca,'XTickLabels',[],'visible','off'); 
annotation('textbox','string','T_{1}^{f}','FontSize',25,'FontWeight','bold','Color','w','LineStyle','none','Position',[0.1,0.8,0.040624998944501,0.06954688970608]);

% subplot(2,2,2); imagesc(squeeze(K_map(:,:,Slice))); axis image;
% colormap('magma'); cb = colorbar; cb.FontSize = 20; cb.Color = 'w'; title(cb,'s^{-1}','FontSize',22,'Color','w'); caxis([0 100]);
% set(gca,'YTickLabels',[],'visible','off'); set(gca,'XTickLabels',[],'visible','off'); 
% annotation('textbox','string','\it{k}','FontSize',28,'FontWeight','bold','Color','w','LineStyle','none');

subplot(3,3,5); imagesc(imgaussfilt(squeeze(M0s_map(67,:,:)),Sigma)); axis image; 
colormap('magma'); cb = colorbar; cb.FontSize = 18; cb.Color = 'w'; caxis([0 0.09]);
set(gca,'YTickLabels',[],'visible','off'); set(gca,'XTickLabels',[],'visible','off');
annotation('textbox','string','\it{f}','FontSize',30,'FontWeight','bold','Color','w','LineStyle','none','Position',[0.1,0.5,0.036458332402011,0.061116963719166]);

subplot(3,3,8); imagesc(imgaussfilt(squeeze(delta_map(67,:,:)),Sigma)); axis image;
colormap('magma'); cb = colorbar; cb.FontSize = 18; cb.Color = 'w'; caxis([0 0.0062]); title(cb,'s','FontSize',20,'Color','w'); 
set(gca,'YTickLabels',[],'visible','off'); set(gca,'XTickLabels',[],'visible','off'); 
annotation('textbox','string','T_{1D}^{s}','FontSize',25,'FontWeight','bold','Color','w','LineStyle','none','Position',[0.1,0.2,0.050520831982916,0.06954688970608]);

subplot(3,3,3); imagesc(imgaussfilt(squeeze(T1_map(:,73,:)),Sigma)); axis image; 
colormap('magma'); cb = colorbar; cb.FontSize = 18; cb.Color = 'w'; title(cb,'s','FontSize',20,'Color','w'); caxis([0 2.5]);
set(gca,'YTickLabels',[],'visible','off'); set(gca,'XTickLabels',[],'visible','off');

subplot(3,3,6); imagesc(imgaussfilt(squeeze(M0s_map(:,73,:)),Sigma)); axis image; 
colormap('magma'); cb = colorbar; cb.FontSize = 18; cb.Color = 'w'; caxis([0 0.09]);
set(gca,'YTickLabels',[],'visible','off'); set(gca,'XTickLabels',[],'visible','off');

subplot(3,3,9); imagesc(imgaussfilt(squeeze(delta_map(:,73,:)),Sigma)); axis image; 
colormap('magma'); cb = colorbar; cb.FontSize = 18; cb.Color = 'w'; caxis([0 0.0062]); title(cb,'s','FontSize',20,'Color','w'); 
set(gca,'YTickLabels',[],'visible','off'); set(gca,'XTickLabels',[],'visible','off');

subplot(3,3,1); imagesc(imgaussfilt(squeeze(T1_map(:,:,37)),Sigma)); axis image; 
colormap('magma'); cb = colorbar; cb.FontSize = 18; cb.Color = 'w'; title(cb,'s','FontSize',20,'Color','w'); caxis([0 2.5]);
set(gca,'YTickLabels',[],'visible','off'); set(gca,'XTickLabels',[],'visible','off');

subplot(3,3,4); imagesc(imgaussfilt(squeeze(M0s_map(:,:,37)),Sigma)); axis image; 
colormap('magma'); cb = colorbar; cb.FontSize = 18; cb.Color = 'w'; caxis([0 0.09]);
set(gca,'YTickLabels',[],'visible','off'); set(gca,'XTickLabels',[],'visible','off');

subplot(3,3,7); imagesc(imgaussfilt(squeeze(delta_map(:,:,37)),Sigma)); axis image; 
colormap('magma'); cb = colorbar; cb.FontSize = 18; cb.Color = 'w'; caxis([0 0.0062]); title(cb,'s','FontSize',20,'Color','w'); 
set(gca,'YTickLabels',[],'visible','off'); set(gca,'XTickLabels',[],'visible','off');

set(gcf,'color','k','InvertHardCopy','off'); colormap(magma());

%% Single voxel signal comparison. Commented out sections are for phantom data.

Slice = 67; %Slice = 28; 
subplot(1,2,1); imagesc(squeeze(M0s_map(Slice,:,:))); axis image; 
colormap('magma'); cb = colorbar; cb.FontSize = 16; cb.Color = 'w'; caxis([0 0.2]);
set(gca,'YTickLabels',[],'visible','off'); set(gca,'XTickLabels',[],'visible','off');
subplot(1,2,2); imagesc(squeeze(delta_map(Slice,:,:))); axis image; 
colormap('magma'); cb = colorbar; cb.FontSize = 16; cb.Color = 'w'; caxis([0 0.0065]);
set(gca,'YTickLabels',[],'visible','off'); set(gca,'XTickLabels',[],'visible','off');

% WM pixel.
%WM_params = squeeze(params_reshaped(15,27,Slice,:));
%WM_data = squeeze(abs(MRF_LRI5(15,27,Slice,:)));
WM_params = squeeze(params_reshaped(Slice,74,71,:));
WM_data = squeeze(abs(MRF_LRI5(Slice,74,71,:)));

% GM pixel.
%GM_params = squeeze(params_reshaped(17,49,Slice,:));
%GM_data = squeeze(abs(MRF_LRI5(17,49,Slice,:)));
GM_params = squeeze(params_reshaped(Slice,96,18,:)); %84,100
GM_data = squeeze(abs(MRF_LRI5(Slice,96,18,:)));

% CSF pixel.
%CSF_params = squeeze(params_reshaped(13,11,Slice,:));
%CSF_data = squeeze(abs(MRF_LRI5(13,11,Slice,:)));
CSF_params = squeeze(params_reshaped(Slice,58,50,:));
CSF_data = squeeze(abs(MRF_LRI5(Slice,58,50,:)));

% Sequence set-up.
flips = d2r(29.51); TR = 5.33e-3; Dur = 2.51e-3; Delta = 8058.48; n1B = 300; nMB = 300;
TBW = 2; B1rms = 4;
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
[G_3B_WM,wloc_3B_WM] = SuperLorentzian_LSint(WM_params(9),df_3B);
[G_2B_WM,wloc_2B_WM] = SuperLorentzian_LSint(WM_params(9),df_2B);
[G_1B_WM,wloc_1B_WM] = SuperLorentzian_LSint(WM_params(9),df_1B);
[G_3B_GM,wloc_3B_GM] = SuperLorentzian_LSint(GM_params(9),df_3B);
[G_2B_GM,wloc_2B_GM] = SuperLorentzian_LSint(GM_params(9),df_2B);
[G_1B_GM,wloc_1B_GM] = SuperLorentzian_LSint(GM_params(9),df_1B);
[G_3B_CSF,wloc_3B_CSF] = SuperLorentzian_LSint(CSF_params(9),df_3B);
[G_2B_CSF,wloc_2B_CSF] = SuperLorentzian_LSint(CSF_params(9),df_2B);
[G_1B_CSF,wloc_1B_CSF] = SuperLorentzian_LSint(CSF_params(9),df_1B);

% Simulate signals.
GM_M0 = [GM_params(1), GM_params(2).*(1-GM_params(3)), GM_params(2).*GM_params(3)]; 
Mss_GM = Dictionary_function_PR_v2(flips,TR,Dur,GM_params(10),GM_params(5:7)',GM_params(8:9)',GM_M0(:),GM_params(4),Delta,TBW,nMB,n1B,b1sqrd,G_1B_GM,G_2B_GM,G_3B_GM,wloc_1B_GM,wloc_2B_GM,wloc_3B_GM);
WM_M0 = [WM_params(1), WM_params(2).*(1-WM_params(3)), WM_params(2).*WM_params(3)]; 
Mss_WM = Dictionary_function_PR_v2(flips,TR,Dur,WM_params(10),WM_params(5:7)',WM_params(8:9)',WM_M0(:),WM_params(4),Delta,TBW,nMB,n1B,b1sqrd,G_1B_WM,G_2B_WM,G_3B_WM,wloc_1B_WM,wloc_2B_WM,wloc_3B_WM);
CSF_M0 = [CSF_params(1), CSF_params(2).*(1-CSF_params(3)), CSF_params(2).*CSF_params(3)]; 
Mss_CSF = Dictionary_function_PR_v2(flips,TR,Dur,CSF_params(10),CSF_params(5:7)',CSF_params(8:9)',CSF_M0(:),CSF_params(4),Delta,TBW,nMB,n1B,b1sqrd,G_1B_CSF,G_2B_CSF,G_3B_CSF,wloc_1B_CSF,wloc_2B_CSF,wloc_3B_CSF);

cm = lines(3);
figure(2);
plot(abs(Mss_WM)./mean(abs(Mss_WM)),'-','LineWidth',3,'Color',cm(1,:)); hold on;
plot(abs(Mss_GM)./mean(abs(Mss_GM)),'-','LineWidth',3,'Color',cm(2,:))
plot(abs(Mss_CSF)./mean(abs(Mss_CSF)),'-','LineWidth',3,'Color',cm(3,:))
plot(1:20:1200,WM_data(1:20:1200)./mean(WM_data(1:20:1200)),'o','LineWidth',3,'Color',cm(1,:));
plot(1:20:1200,GM_data(1:20:1200)./mean(GM_data(1:20:1200)),'o','LineWidth',3,'Color',cm(2,:));
plot(1:20:1200,CSF_data(1:20:1200)./mean(CSF_data(1:20:1200)),'o','LineWidth',3,'Color',cm(3,:));
ll = legend('WM','GM','CSF'); legend boxoff; ll.FontSize = 16; ll.Orientation = 'horizontal';
set(gca,'FontSize',16); grid on; grid minor;
xlabel('Pulse No.','FontSize',16); ylabel('Norm. Signal (a.u.)','FontSize',16);
set(gcf, 'Position',  [100, 100, 1467, 528])