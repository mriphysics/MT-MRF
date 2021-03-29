%% Calculation of signal differentials and CRLB. Daniel West 2020.

close all; clear all; clc

plot_case = 'PNRs'; % Decide whether to plot Fig3 ('PNRs') or signal differentials (not included in manuscript).

% Define default tissue parameters.
T1f = 1.36; T1s = 1; T1d = 6.05; T2f = 84; T2s = 8.28e-3; % First two in s and latter three in ms.
delta = 1; k = 55.2; M0 = 1; M0s = 0.0679;
x0 = [T1f,T1s,T1d,T2f,T2s,M0s,M0,k,delta];

% Define acquisition settings.
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

df_1B = [0 0 0]; df_2B = [0 0 Delta]; df_3B = [-Delta 0 Delta];

[G_3B,wloc_3B] = SuperLorentzian_LSint(T2s*1e-3,df_3B);
[G_2B,wloc_2B] = SuperLorentzian_LSint(T2s*1e-3,df_2B);
[G_1B,wloc_1B] = SuperLorentzian_LSint(T2s*1e-3,df_1B);

% Calculate default signal for scaling.
sig_default = Dictionary_function_CRLB(x0,flips,TR,Dur,dphi,Delta,TBW,nMB,n1B,b1sqrd,G_1B,G_2B,G_3B,wloc_1B,wloc_2B,wloc_3B);

%% Generate parameter combinations.

Npar = length(x0); Nmeas = (2*(n1B + nMB));

Sigma = 0.004;
mean_sig = mean(abs(sig_default));
SNR = mean_sig/Sigma;

% Different parameter combinations.
Param_List = 1:9;

Param_List_9 = nchoosek(Param_List,9);
Param_List_8 = nchoosek(Param_List,8);
Param_List_7 = nchoosek(Param_List,7);
Param_List_6 = nchoosek(Param_List,6);
Param_List_5 = nchoosek(Param_List,5);
Param_List_4 = nchoosek(Param_List,4);
Param_List_3 = nchoosek(Param_List,3);
Param_List_2 = nchoosek(Param_List,2);
Param_List_1 = nchoosek(Param_List,1);

PL_comp = zeros(9+36+84+126+126+84+36+9+1,9);
PL_comp(1:9,1) = Param_List_1;
PL_comp(10:10+35,1:2) = Param_List_2;
PL_comp(46:46+83,1:3) = Param_List_3;
PL_comp(130:130+125,1:4) = Param_List_4;
PL_comp(256:256+125,1:5) = Param_List_5;
PL_comp(382:382+83,1:6) = Param_List_6;
PL_comp(466:466+35,1:7) = Param_List_7;
PL_comp(502:502+8,1:8) = Param_List_8;
PL_comp(511,:) = Param_List_9;

switch plot_case
    case 'diffs'
        PL_comp = Param_List_9;
    case 'PNRs'
end

%% Calculate differentials and CRLB.

PNR = zeros(size(PL_comp)); GT = zeros(size(PL_comp));

for tt = 1:size(PL_comp,1)
    
    dfdx = zeros(Npar,Nmeas); delta_s = zeros(Npar,Nmeas); dfdx_rel = zeros(Npar,Nmeas);
    nz_elem = nonzeros(PL_comp(tt,:)); nz_size = size(nz_elem,1);

    for pp = nz_elem'
        
        x0_plus = x0; x0_minus = x0;
        
        dx = 0.005*x0(pp);
        
        x0_plus(pp) = x0_plus(pp) + dx;
        x0_minus(pp) = x0_minus(pp) - dx;
        
        if pp == 5
            
            % Recompute saturation terms.
            [G_3B_p,wloc_3B_p] = SuperLorentzian_LSint(x0_plus(pp)*1e-3,df_3B);
            [G_2B_p,wloc_2B_p] = SuperLorentzian_LSint(x0_plus(pp)*1e-3,df_2B);
            [G_1B_p,wloc_1B_p] = SuperLorentzian_LSint(x0_plus(pp)*1e-3,df_1B);
            
            [G_3B_m,wloc_3B_m] = SuperLorentzian_LSint(x0_minus(pp)*1e-3,df_3B);
            [G_2B_m,wloc_2B_m] = SuperLorentzian_LSint(x0_minus(pp)*1e-3,df_2B);
            [G_1B_m,wloc_1B_m] = SuperLorentzian_LSint(x0_minus(pp)*1e-3,df_1B);
            
            sig_plus = Dictionary_function_CRLB(x0_plus,flips,TR,Dur,dphi,Delta,TBW,nMB,n1B,b1sqrd,G_1B_p,G_2B_p,G_3B_p,wloc_1B_p,wloc_2B_p,wloc_3B_p);
            sig_minus = Dictionary_function_CRLB(x0_minus,flips,TR,Dur,dphi,Delta,TBW,nMB,n1B,b1sqrd,G_1B_m,G_2B_m,G_3B_m,wloc_1B_m,wloc_2B_m,wloc_3B_m);
            
        else
            
            sig_plus = Dictionary_function_CRLB(x0_plus,flips,TR,Dur,dphi,Delta,TBW,nMB,n1B,b1sqrd,G_1B,G_2B,G_3B,wloc_1B,wloc_2B,wloc_3B);
            sig_minus = Dictionary_function_CRLB(x0_minus,flips,TR,Dur,dphi,Delta,TBW,nMB,n1B,b1sqrd,G_1B,G_2B,G_3B,wloc_1B,wloc_2B,wloc_3B);
            
        end
        
        dfdx(pp,:) = (abs(sig_plus)-abs(sig_minus))/(2*dx);
        delta_s(pp,:) = dfdx(pp,:)*0.05*x0(pp);
        dfdx_rel(pp,:) = (dfdx(pp,:)*x0(pp))./abs(sig_default);
        
    end
    
    [CRLB,FIM] = CRLB_calc(dfdx(nonzeros(PL_comp(tt,:)),:)); CRLB_noise = sqrt(CRLB)*Sigma;
    GT(tt,nz_elem) = x0(nonzeros(PL_comp(tt,:))); PNR(tt,nz_elem) = (GT(tt,nz_elem)'./CRLB_noise)/SNR;
    
end

% Remove unwanted combinations.
ind6 = find(PNR(:,6) == 0);
ind7 = find(PNR(:,7) == 0);
ind3 = find(PNR(:,3) == 0);
ind_673 = unique([ind6;ind7;ind3],'rows');
PNR_matrix_673 = PNR; PNR_matrix_673(ind_673,:) = [];

%% Plot signal differentials for Fig 8.

switch plot_case
    
    case 'diffs'

        title_str = {'T_{1}^{f}','T_{1Z}^{s}','T_{1D}^{s}','T_{2}^{f}','T_{2}^{s}','\it{f}','\it{k}','\delta'};

        figure(1);
        subplot(2,4,1);plot(dfdx_rel(1,:),'-','LineWidth',3,'Color','k');
        set(gca,'FontSize',14); grid on; grid minor; xlabel('Pulse No.','FontSize',14); ylabel('S_{\theta}','FontSize',24);
        dummyh = line(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
        legend(dummyh, 'String',title_str{1},'FontSize',20,'Location','southeast'); legend boxoff;

        subplot(2,4,2);plot(dfdx_rel(2,:),'-','LineWidth',3,'Color','k');
        set(gca,'FontSize',14); grid on; grid minor; xlabel('Pulse No.','FontSize',14); ylabel('S_{\theta}','FontSize',24);
        dummyh = line(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
        legend(dummyh, 'String',title_str{2},'FontSize',20,'Location','southwest'); legend boxoff;

        subplot(2,4,3);plot(dfdx_rel(3,:),'-','LineWidth',3,'Color','k');
        set(gca,'FontSize',14); grid on; grid minor; xlabel('Pulse No.','FontSize',14); ylabel('S_{\theta}','FontSize',24);
        dummyh = line(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
        legend(dummyh, 'String',title_str{3},'FontSize',20); legend boxoff;

        subplot(2,4,4);plot(dfdx_rel(4,:),'-','LineWidth',3,'Color','k');
        set(gca,'FontSize',14); grid on; grid minor; xlabel('Pulse No.','FontSize',14); ylabel('S_{\theta}','FontSize',24);
        dummyh = line(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
        legend(dummyh, 'String',title_str{4},'FontSize',20,'Location','southeast'); legend boxoff;

        subplot(2,4,5);plot(dfdx_rel(5,:),'-','LineWidth',3,'Color','k');
        set(gca,'FontSize',14); grid on; grid minor; xlabel('Pulse No.','FontSize',14); ylabel('S_{\theta}','FontSize',24);
        dummyh = line(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
        legend(dummyh, 'String',title_str{5},'FontSize',20); legend boxoff;

        subplot(2,4,6);plot(dfdx_rel(6,:),'-','LineWidth',3,'Color','k');
        set(gca,'FontSize',14); grid on; grid minor; xlabel('Pulse No.','FontSize',14); ylabel('S_{\theta}','FontSize',24);
        dummyh = line(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
        legend(dummyh, 'String',title_str{6},'FontSize',20,'Location','southeast'); legend boxoff;

        % Differential wrt M0 gives is flat.
        %subplot(3,3,7);plot(dfdx_rel(7,:),'-','LineWidth',3,'Color','k');
        %set(gca,'FontSize',14); grid on; grid minor; xlabel('Pulse No.','FontSize',14); ylabel('S_{\theta}','FontSize',24);
        %dummyh = line(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none'); ylim([0.9 1.1])
        %legend(dummyh, 'String',title_str{7},'FontSize',20); legend boxoff;

        subplot(2,4,7);plot(dfdx_rel(8,:),'-','LineWidth',3,'Color','k');
        set(gca,'FontSize',14); grid on; grid minor; xlabel('Pulse No.','FontSize',14); ylabel('S_{\theta}','FontSize',24);
        dummyh = line(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
        legend(dummyh, 'String',title_str{7},'FontSize',20,'Location','southeast'); legend boxoff;

        subplot(2,4,8);plot(dfdx_rel(9,:),'-','LineWidth',3,'Color','k');
        set(gca,'FontSize',14); grid on; grid minor; xlabel('Pulse No.','FontSize',14); ylabel('S_{\theta}','FontSize',24);
        dummyh = line(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
        legend(dummyh, 'String',title_str{8},'FontSize',26); legend boxoff;
        
    case 'PNRs'
        
        PNR_matrix_plot = PNR_matrix_673(:,[6,3,1,2,4,5,8,9]).';
        PNR_matrix_plot = [PNR_matrix_plot;zeros(1,64)];
        PNR_matrix_plot = [PNR_matrix_plot,zeros(9,1)];
        
        figure(2)
        ax1 = axes;
        pcolor(ax1,PNR_matrix_plot); set(gca,'FontSize',22); caxis([0 2]);
        xlabel('Combination No.','Position',[32.2,0.55,0]); 
        set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
        annotation('textbox','String','\delta','FontSize',30,'LineStyle','none','Position',[0.06,0.133,0.0333,0.057]);
        annotation('textbox','String','\it{k}','FontSize',30,'LineStyle','none','Position',[0.06,0.233,0.0333,0.057]);
        annotation('textbox','String','T_{2}^{s}','FontSize',24,'LineStyle','none','Position',[0.06,0.353,0.0333,0.057]);
        annotation('textbox','String','T_{2}^{f}','FontSize',24,'LineStyle','none','Position',[0.06,0.453,0.0333,0.057]);
        annotation('textbox','String','T_{1Z}^{s}','FontSize',24,'LineStyle','none','Position',[0.06,0.553,0.0333,0.057]);
        annotation('textbox','String','T_{1}^{f}','FontSize',24,'LineStyle','none','Position',[0.06,0.653,0.0333,0.057]);
        annotation('textbox','String','T_{1D}^{s}','FontSize',24,'LineStyle','none','Position',[0.06,0.763,0.0333,0.057]);
        annotation('textbox','String','\it{f}','FontSize',30,'LineStyle','none','Position',[0.06,0.843,0.0333,0.057]);
        annotation('textbox','String','10','FontSize',24,'LineStyle','none','Position',[0.1980,0.056,0.0333,0.057]);
        annotation('textbox','String','20','FontSize',24,'LineStyle','none','Position',[0.3180,0.056,0.0333,0.057]);
        annotation('textbox','String','30','FontSize',24,'LineStyle','none','Position',[0.4360,0.056,0.0333,0.057]);
        annotation('textbox','String','40','FontSize',24,'LineStyle','none','Position',[0.5540,0.056,0.0333,0.057]);
        annotation('textbox','String','50','FontSize',24,'LineStyle','none','Position',[0.6720,0.056,0.0333,0.057]);
        annotation('textbox','String','60','FontSize',24,'LineStyle','none','Position',[0.790,0.056,0.0333,0.057]);
        ax2 = axes;
        pcolor(ax2,PNR_matrix_plot); set(gca,'FontSize',22); caxis([0 2]); set(gca,'YDir','reverse');
        linkaxes([ax1,ax2])
        ax2.Visible = 'off'; ax2.XTick = []; ax2.YTick = [];
        colormap(ax2,[[1 1 1];viridis(10000)]);
        colormap(ax1,viridis(10000));
        set([ax1,ax2],'Position',[.1 .11 .755 .815]);
        cb2 = colorbar(ax1,'Position',[.88 .11 .0275 .815]);
        title(cb2,'PNR/SNR','FontSize',24); 
        cb2.Ticks = [0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0]; cb2.TickLabels = {'0','0.25', '0.5','0.75', '1.0','1.25','1.5','1.75',''};

end
