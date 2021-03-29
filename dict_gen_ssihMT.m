%%% ss-ihMT dictionary generation. DW16 2021 script for phantom validation. %%%

%%% Relies on scripts from: https://github.com/mriphysics/ihMT_steadystate .

close all; clear all; clc;

%% Generate parameter combinations.

n_increments = 100;
T1f = linspace(0.2,4,n_increments);
delta = l1;
M0s = linspace(0,0.3,n_increments);
T1s = 1;
T1d = linspace(0,40e-3,n_increments);
T2f = 84e-3; %130e-3;
T2s = 8.28; % This is now MICROSECONDS.
k = 55.2;
B1_sf = 1;

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
Var_Params_noMT = unique(combvec(M0_noMT.',T12x.',B1_sf).','rows');
Var_Params_MT = unique(combvec(M0_MT.',T12x.',B1_sf).','rows');
Var_Params_ihMT = unique(combvec(M0_ihMT.',T12x.',B1_sf).','rows');

% Sort to order: [R1f R2 M0s  R1sz  R1D  f  T2s(us)  k   sf   b1sf = 1] -->
% CHECK SUBSEQUENT FUNCTIONS - 3RD PARAM IS ACTUALLY M0s (NOT M0f).
param_matrix = [Var_Params_noMT;Var_Params_MT;Var_Params_ihMT];
clear Var_Params_noMT Var_Params_MT Var_Params_ihMT
param_matrix_mod = param_matrix(:,[5 8 2 6 7 3 9 4 1 10]);

% Convert relaxation times into rates.
param_matrix_mod(:,[1 2 4 5]) = 1./param_matrix_mod(:,[1 2 4 5]);

% This creates some Inf, so assign a very large number to these values.
[idx,idy] = find(isinf(param_matrix_mod));
param_matrix_mod(idx,idy) = 100000;

% Ninth column was M0s but now should be a scaling factor.
param_matrix_mod(:,9) = -1; % Forces normalisation in fwd_model functions.

%% Define SSFP pulse sequence (SJM).

seq_pars = struct;
seq_pars.flips = deg2rad(10:10:80);
nf = length(seq_pars.flips);
seq_pars.tau = 2.2e-3;
seq_pars.TR = 5e-3;
seq_pars.b1_rms = 4.15;
seq_pars.delta = 8e3;
dt = 10e-6;

seq_pars.b1sqrd = {};

for ii=1:nf
    % 3 bands
    [~,seq_pars.b1sqrd{ii,3},~,seq_pars.TBP] = gen_MB_pulse(seq_pars.flips(ii),seq_pars.tau,seq_pars.TR,seq_pars.b1_rms,seq_pars.delta,'3','alpha',3,'dt',dt); 
    % 2 bands
    [~,seq_pars.b1sqrd{ii,2},~] = gen_MB_pulse(seq_pars.flips(ii),seq_pars.tau,seq_pars.TR,seq_pars.b1_rms,seq_pars.delta,'2+','alpha',3,'dt',dt);
    % 1 bands
    seq_pars.b1sqrd{ii,1} = seq_pars.b1sqrd{ii,2}.*[0 1 0]; 
end
seq_pars.Delta_Hz = seq_pars.delta * [-1 0 1];
seq_pars.dphi=0;

%% Define SPGR pulse sequence (SJM).

seq_pars_spgr = seq_pars;
seq_pars_spgr.flips = deg2rad(2:2:16);
nf = length(seq_pars_spgr.flips);

seq_pars_spgr.b1sqrd = {};

for ii=1:nf
    % 3 bands
    [~,seq_pars_spgr.b1sqrd{ii,3},~,TBP] = gen_MB_pulse(seq_pars_spgr.flips(ii),seq_pars_spgr.tau,seq_pars_spgr.TR,seq_pars_spgr.b1_rms,seq_pars_spgr.delta,'3','alpha',3,'dt',dt); 
    % 2 bands
    [~,seq_pars_spgr.b1sqrd{ii,2},~] = gen_MB_pulse(seq_pars_spgr.flips(ii),seq_pars_spgr.tau,seq_pars_spgr.TR,seq_pars_spgr.b1_rms,seq_pars_spgr.delta,'2+','alpha',3,'dt',dt);
    % 1 bands
    seq_pars_spgr.b1sqrd{ii,1} = seq_pars_spgr.b1sqrd{ii,2}.*[0 1 0]; 
end
seq_pars_spgr.Delta_Hz = seq_pars.delta * [-1 0 1];

% Add echo time.
seq_pars_spgr.TE = 2.5e-3;

%% Generate signals for dictionary.

[G,wloc] = SuperLorentzian_LSint(T2s*1e-6,seq_pars.Delta_Hz);

flag = []; % Empty for super-Lorentzian.

delete(gcp('nocreate')); c = parcluster('local'); c.NumWorkers = 16; parpool(c, c.NumWorkers);

ssfp_dict = zeros(size(param_matrix_mod,1),8,3);
spgr_dict = zeros(size(param_matrix_mod,1),8,3);
tic
parfor jj = 1:size(param_matrix_mod,1)
    ssfp_dict(jj,:,:) = ssfp_ihmt_fit_fwd_model(param_matrix_mod(jj,:),seq_pars,flag,G,wloc);
    spgr_dict(jj,:,:) = spgr_ihmt_fit_fwd_model(param_matrix_mod(jj,:),seq_pars_spgr,flag,G,wloc);
end
toc
dict_all = cat(3,ssfp_dict,spgr_dict); clear ssfp_dict spgr_dict

% Combine second and third dimensions.
dict_all_rs = reshape(dict_all,[size(dict_all,1) size(dict_all,2)*size(dict_all,3)]); clear dict_all;

clearvars -except dict_all_rs param_matrix_mod

%% Load in data to fit.

load bin/fitdata.mat
IX = 3; % PL161
IX = 2; % BSA
IX = 5; % MnCl2

IX = [5 2 3];

% Get the data for each phantom into a cell array.
xdata = {};
sdata = {};
for jj=1:3
    
    tmp = squeeze(kdata{IX(jj),1});
    % Average the 2+ and 2- data.
    tmp(:,:,2) = 0.5*(tmp(:,:,2)+tmp(:,:,3));
    tmp(:,:,3) = [];
    
    % Now SPGR data.
    tmp2 = squeeze(kdata{IX(jj),2});
    tmp2(:,:,2) = 0.5*(tmp2(:,:,2)+tmp2(:,:,3));
    tmp2(:,:,3) = [];
    
    % Reorder and scale down by 1000. Change?
    xdata{jj,1} = permute(tmp/1000,[2 3 1]);
    xdata{jj,2} = permute(tmp2/1000,[2 3 1]);
 
    % Standard deviation.
    sdata{jj,1} = std(xdata{jj,1},1,3);
    sdata{jj,2} = std(xdata{jj,2},1,3);
    % Store number of samples.
    sdata{jj,3} = size(xdata{jj,1},3);
    
end

% Separate-out xdata.
MnCl_ssfp = xdata{1,1}; BSA_ssfp = xdata{2,1}; PL161_ssfp = xdata{3,1};
MnCl_spgr = xdata{1,2}; BSA_spgr = xdata{2,2}; PL161_spgr = xdata{3,2};

% Combine sequences.
MnCl_data = cat(2,MnCl_ssfp,MnCl_spgr); BSA_data = cat(2,BSA_ssfp,BSA_spgr); PL161_data = cat(2,PL161_ssfp,PL161_spgr);

% Combine phantoms.
alldata_phantom = cat(3,MnCl_data,BSA_data,PL161_data);

% Permute dimensions (pixels x FAs x datasets).
alldata_mod = permute(alldata_phantom,[3 1 2]);

% Combine second and third dimensions.
alldata_mod_rs = reshape(alldata_mod,[size(alldata_mod,1) size(alldata_mod,2)*size(alldata_mod,3)]);
alldata_mod_rs = permute(alldata_mod_rs, [2 1]);

%% Perform fitting.

% Normalise dictionary again?
norm_atoms = conj(dict_all_rs./sum(abs(dict_all_rs).^2,2).^0.5);

sig_norm = conj(alldata_mod_rs./sum(abs(alldata_mod_rs).^2,1).^0.5);

nchunks = ceil(size(norm_atoms,1)*size(sig_norm,2)/(30*1024*1024*128));
npos = size(alldata_mod_rs,2);
 
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
param_fit = param_matrix_mod(max_idx(:),:);

%% Analyse results.

MnCl2_estimates = param_fit(1:size(MnCl_data,3),:);
BSA_estimates = param_fit(size(MnCl_data,3)+1:size(MnCl_data,3)+size(BSA_data,3),:);
PL161_estimates = param_fit(size(MnCl_data,3)+size(BSA_data,3)+1:end,:);

% Extract parameters-of-interest.
ss_T1f_MnCl2 = 1./MnCl2_estimates(:,1);
ss_T1f_BSA = 1./BSA_estimates(:,1);
ss_T1f_PL161 = 1./PL161_estimates(:,1);

ss_f_MnCl2 = MnCl2_estimates(:,3);
ss_f_BSA = BSA_estimates(:,3);
ss_f_PL161 = PL161_estimates(:,3);

ss_T1D_MnCl2 = 1./MnCl2_estimates(:,5);
ss_T1D_BSA = 1./BSA_estimates(:,5);
ss_T1D_PL161 = 1./PL161_estimates(:,5);

% Key values.
f_BSA = mean(ss_f_BSA); fsd_BSA = std(ss_f_BSA);
f_PL161 = mean(ss_f_PL161); fsd_PL161 = std(ss_f_PL161);
dipolar_PL161 = mean(ss_T1D_PL161); dipolarsd_PL161 = std(ss_T1D_PL161);
