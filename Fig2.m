%% Generates efficiency heat maps for Figure 2. Daniel West 2020

close all; clear all; clc

% 1.5T actual parameter set.
T1x = [650 1000 6.5]*1e-3; T2x = [80e-3, 12.5e-6]; delta = 0.65;
k = 65; M0s_Varma = (7.3*(1/T1x(1)))/k; M0s = M0s_Varma/(1 + M0s_Varma); 
M0f = 1-M0s; M0 = [M0f M0s*(1-delta) M0s*delta];
ihMT_Params = [29.51,5.33,2.51,8.05848,300,300]; dphi = 0; TBP = 2; B1rms = 4.0;

% Evaluate signal function for different df/FA combinations.
FA_var = linspace(1,50,50); df_var = linspace(3,12,50);

Eff_FAdf = zeros(length(FA_var),length(df_var)); Con_FAdf = zeros(length(FA_var),length(df_var));
for ii = 1:length(FA_var)
    disp(['FA = ',num2str(FA_var(ii))])
    for jj = 1:length(df_var)
        Eff_FAdf(ii,jj) = ssSSFP_ihMT_Pulses_2BP(FA_var(ii),ihMT_Params(2),ihMT_Params(3),dphi,T1x,T2x,M0,k,df_var(jj),TBP,ihMT_Params(5),ihMT_Params(6),B1rms,'MRF');
        Con_FAdf(ii,jj) = abs(Eff_FAdf(ii,jj))*sqrt(ihMT_Params(2)*1e-3);
    end
end

% Evaluate signal function for different timing combinations.
Dur_var = linspace(1.5,4,50); TR_var = linspace(4.5,8,50);

Eff_DurTR = zeros(length(Dur_var),length(TR_var)); Con_DurTR = zeros(length(Dur_var),length(TR_var));
for ii = 1:length(Dur_var)
    disp(['Pulse Duration = ',num2str(Dur_var(ii))])
    for jj = 1:length(TR_var)
        Eff_DurTR(ii,jj) = ssSSFP_ihMT_Pulses_2BP(ihMT_Params(1),TR_var(jj),Dur_var(ii),dphi,T1x,T2x,M0,k,ihMT_Params(4),TBP,ihMT_Params(5),ihMT_Params(6),B1rms,'MRF');
        Con_DurTR(ii,jj) = abs(Eff_DurTR(ii,jj))*sqrt(TR_var(jj)*1e-3); 
    end
end

% Evaluate signal function for different pulse number combinations.
nMB_var = 2:2:400; nSB_var = 2:2:400; B1max_CalcTest = zeros(length(nMB_var),length(nSB_var));

Eff_npulses = zeros(length(nMB_var),length(nSB_var)); Con_npulses = zeros(length(nMB_var),length(nSB_var));
for ii = 1:length(nMB_var)
    disp(['no. MB = ',num2str(nMB_var(ii))])
    for jj = 1:length(nSB_var)
        Eff_npulses(ii,jj) = ssSSFP_ihMT_Pulses_2BP(ihMT_Params(1),ihMT_Params(2),ihMT_Params(3),dphi,T1x,T2x,M0,k,ihMT_Params(4),TBP,nMB_var(ii),nSB_var(jj),B1rms,'MRF');
        [~,~,~,~,B1max_CalcTest(ii,jj)] = gen_CSMT_pulse_Diffnp(d2r(ihMT_Params(1)),ihMT_Params(3)*0.001,ihMT_Params(2)*0.001,B1rms,ihMT_Params(4)*1000,3,nMB_var(ii),nSB_var(jj),'sigma',2);
        Con_npulses(ii,jj) = abs(Eff_npulses(ii,jj))*sqrt(ihMT_Params(2)*1e-3);    
    end
end

%% Plot results.

ylims = [0.11 0.815]; xlims = [0.08 0.2]; dx = 0.3;

colormap(gray)
figure(1);
set(gcf,'Units','normalized','Outerposition',[0 0 1 0.8],'Color','w');
subplot(1,3,1);
imagesc(df_var,FA_var,abs(Eff_FAdf)); axis square; hold on;
plot(ihMT_Params(4),ihMT_Params(1),'.','Color',[0.9290    0.6940    0.1250],'MarkerSize',40)
set(gca,'FontSize',18); xlabel('\Delta (kHz)','FontSize',18); ylabel('FA (^{o})','FontSize',18); cb = colorbar; title(cb,'\eta (s^{-1/2})','FontSize',18); caxis([0.05 0.15]); 
sub1_pos = get(gca,'Position'); 
set(gca,'Position',[xlims(1) ylims(1) xlims(2) ylims(2)],'Units','normalized');
annotation('textbox','String','(a)','FontSize',24,'LineStyle','none','Position',[0.08,0.77,0.05,0.07])

figure(1); subplot(1,3,2);
im = imagesc(nSB_var,nMB_var,abs(Eff_npulses)); hold on; axis square; 
% c is found by: imagesc(nSB_var,nMB_var,B1max_CalcTest); caxis([18 20]).
c = area([2 500 500 500],[2 376 500 2],0); d = 0.2; c.FaceAlpha=d; c.FaceColor = [1 0 0]; c.EdgeColor = [1 0 0]; c.EdgeAlpha = d;
d = plot(ihMT_Params(6),ihMT_Params(5),'.','Color',[0.9290    0.6940    0.1250],'MarkerSize',40);
e = plot(152,141,'.','Color',[0.4940    0.1840    0.5560],'MarkerSize',40);
delta = plot(350,100,'.','Color',[0.4660    0.6740    0.1880],'MarkerSize',40);
set(gca,'FontSize',18); xlabel('n_{1B}','FontSize',18); ylabel('n_{MB}','FontSize',18); cb = colorbar; title(cb,'\eta (s^{-1/2})','FontSize',18); caxis([0.05 0.16])
sub2_pos = get(gca,'Position');
ll = legend([d,e,delta,c], {'Scan','Optimum','Forbidden','B_{1,max}'},'Location','northoutside'); 
ll.FontSize = 18; legend boxoff; ll.Orientation = 'vertical';
set(gca,'Position',[xlims(1)+dx ylims(1) xlims(2) ylims(2)],'Units','normalized');
annotation('textbox','String','(b)','FontSize',24,'LineStyle','none','Position',[0.38,0.77,0.05,0.07])
xticks([100 200 300 400]); yticks([100 200 300 400]); xticklabels({'100','200','300','400'}); yticklabels({'100','200','300','400'})

subplot(1,3,3);
imagesc(TR_var,Dur_var,abs(Eff_DurTR)); axis square; hold on;
a = area(2.1*Dur_var,Dur_var,10); b = 0.2; a.FaceAlpha=b; a.FaceColor = [0 1 0]; a.EdgeColor = [0 1 0]; a.EdgeAlpha = b;
e = area(2.5+Dur_var,Dur_var,10); delta = 0.2; e.FaceAlpha=delta; e.FaceColor = [0 0 1]; e.EdgeColor = [0 0 1]; e.EdgeAlpha = b;
plot(ihMT_Params(2),ihMT_Params(3),'.','Color',[0.9290    0.6940    0.1250],'MarkerSize',40)
set(gca,'FontSize',18); xlabel('TR (ms)','FontSize',18); ylabel('\tau (ms)','FontSize',18); cb = colorbar; title(cb,'\eta (s^{-1/2})','FontSize',18); caxis([0.12 0.15]); 
sub3_pos = get(gca,'Position');
ll = legend([a,e], 'Duty Cycle','Timing','Location','northoutside'); ll.FontSize = 18; legend boxoff;
set(gca,'Position',[xlims(1)+(2*dx) ylims(1) xlims(2) ylims(2)],'Units','normalized');
annotation('textbox','String','(c)','FontSize',24,'LineStyle','none','Position',[0.68,0.77,0.05,0.07])
