%% Extract from script that handles phantom and in vivo recon. data. Daniel West 2020

close all; clear all; clc

%% Phantom data processing for Figure 5.

close all; clear all; clc;

% This is the latest balanced data.
load('/home/dw16/Desktop/MRF_ihMT_Data/Oblique_Phantom/hdp3d_ext.mat')

% Load dictionary.
load('/home/dw16/Desktop/MRF_ihMT_Data/Oblique_Phantom/Phantom_Dict_ex2.mat'); Ur = V(:,1:5);

% Transform phantom data.
Slice = 50;

MRF_sing = squeeze(hdp3d(:,:,Slice,:));
mrf_aux = reshape(MRF_sing,[size(MRF_sing,1)*size(MRF_sing,2) size(MRF_sing,3)]); mrf_aux = mrf_aux';
MRF_LRI5 = (Ur*mrf_aux)'; MRF_LRI5 = reshape(MRF_LRI5,[size(MRF_sing,1) size(MRF_sing,2) size(MRF_LRI5,2)]);

% Take out global phase term.
phase = repmat(angle(sum(MRF_LRI5,3)),[1 1 1200]);
MRF_LRI5 = MRF_LRI5.*exp(-1i*phase);

% Calculate ratios using multiple time-points.
SB_av = (MRF_LRI5(:,:,600)+MRF_LRI5(:,:,1200))./2;
TwoB_av = MRF_LRI5(:,:,300); ThreeB_av = MRF_LRI5(:,:,900);
ihMT(:,:) = real((TwoB_av-ThreeB_av)./(SB_av));
MT(:,:) = 1-(real((TwoB_av+ThreeB_av)./(2*SB_av)));

% Load balanced 1DFT data.
load('/home/dw16/Desktop/MRF_ihMT_Data/Oblique_Phantom/09012020_1D_Data.mat')
OneDFT_Data = abs(imgs{1,4});
time_course = 3601-7:4800-7;
PL161_1DFT = sum(OneDFT_Data(34:41,time_course),1)./mean(sum(OneDFT_Data(34:41,time_course),1));
BSA_1DFT = sum(OneDFT_Data(53:59,time_course),1)./mean(sum(OneDFT_Data(53:59,time_course),1));
Water_1DFT = sum(OneDFT_Data(19:23,time_course),1)./mean(sum(OneDFT_Data(19:23,time_course),1));

% Format balanced MTF data.
BSA_signal = squeeze(mean(abs(MRF_LRI5(73:77,92:97,:)),[1 2]))./mean(mean(abs(MRF_LRI5(73:77,92:97,:)),[1 2]));
BSA_std = squeeze(std(abs(MRF_LRI5(73:77,92:97,:)),0,[1 2])./mean(mean(abs(MRF_LRI5(73:77,92:97,:)),[1 2])));
PL161_signal = squeeze(mean(abs(MRF_LRI5(72:82,72:81,:)),[1 2]))./mean(mean(abs(MRF_LRI5(72:82,72:81,:)),[1 2]));
PL161_std = squeeze(std(abs(MRF_LRI5(72:82,72:81,:)),0,[1 2])./mean(mean(abs(MRF_LRI5(72:82,72:81,:)),[1 2])));
MnCl2_signal = squeeze(mean(abs(MRF_LRI5(73:77,58:62,:)),[1 2]))./mean(mean(abs(MRF_LRI5(73:77,58:62,:)),[1 2]));
MnCl2_std = squeeze(std(abs(MRF_LRI5(73:77,58:62,:)),0,[1 2])./mean(mean(abs(MRF_LRI5(73:77,58:62,:)),[1 2])));

% Create mask.
mask = mean(abs(MRF_LRI5(45:105,50:110)),3) > 1.5;

figure(1); cm = lines(5);
cm_default = magma;
subplot(1,2,2); imagesc((abs(ihMT(45:105,50:110)).*mask)*100); axis image; set(gca,'YTickLabels',[]); set(gca,'XTickLabels',[]); colormap(gca,[[0 0 0]; magma]); caxis([0 30]); annotation('textbox','String','ihMTR (%)','FontSize',30,'Color','w','LineStyle','none','Position',[0.7,0.8,0.174164957521152,0.109848482139183]); cb = colorbar; cb.FontSize = 18; cb.Color = 'w';  set(gca,'YTickLabels',[],'visible','off'); set(gca,'XTickLabels',[],'visible','off');
subplot(1,2,1); imagesc((abs(MT(45:105,50:110)).*mask)*100); axis image; set(gca,'YTickLabels',[]); set(gca,'XTickLabels',[]); colormap(gca,[[0 0 0]; magma]); caxis([0 70]); annotation('textbox','String','MTR (%)','FontSize',30,'Color','w','LineStyle','none','Position',[0.28,0.8,0.148943417713849,0.109848482139183]); cb = colorbar; cb.Color = 'w'; cb.FontSize = 18;  set(gca,'YTickLabels',[],'visible','off'); set(gca,'XTickLabels',[],'visible','off');
set(gcf, 'Position',  [100, 100, 1920, 528]); 
annotation('textbox','String','(c)','FontSize',32,'LineStyle','none','Color','w','Position',[0.0123,0.8869,0.05,0.07])
annotation('textbox','String','MnCl_{2}-Doped Water','FontSize',20,'LineStyle','none','Color','w','Position',[0.085237968643492,0.350915151515152,0.05,0.07])
annotation('textbox','String','PL161','FontSize',20,'LineStyle','none','Color','w','Position',[0.229068916155419,0.267581818181819,0.05,0.07])
annotation('textbox','String','BSA','FontSize',20,'LineStyle','none','Color','w','Position',[0.32041179277437,0.333869696969698,0.05,0.07])
annotation('textbox','String','MnCl_{2}-Doped Water','FontSize',20,'LineStyle','none','Color','w','Position',[0.525237968643492,0.350915151515152,0.05,0.07])
annotation('textbox','String','PL161','FontSize',20,'LineStyle','none','Color','w','Position',[0.669068916155419,0.267581818181819,0.05,0.07])
annotation('textbox','String','BSA','FontSize',20,'LineStyle','none','Color','w','Position',[0.76041179277437,0.333869696969698,0.05,0.07])
set(gcf,'color','k','InvertHardCopy','off');

figure(2)
subplot(5,2,1); imagesc((abs(hdp3d(65:90,50:110,50,1)))); axis image; set(gca,'YTickLabels',[],'visible','off'); set(gca,'XTickLabels',[],'visible','off'); title('#1','FontSize',20); caxis([0 200]); cb = colorbar; cb.FontSize = 16; cb.Color = 'w';
annotation('textbox','String', 'R = 1','Position',[0.13,0.84,0.07,0.07],'FontSize',28,'Color','w','LineStyle','none','FontWeight','bold');
subplot(5,2,3); imagesc((abs(hdp3d(65:90,50:110,50,2)))); axis image; set(gca,'YTickLabels',[],'visible','off'); set(gca,'XTickLabels',[],'visible','off'); title('#2','FontSize',20); caxis([0 20]); cb = colorbar; cb.FontSize = 16; cb.Color = 'w';
annotation('textbox','String', 'R = 2','Position',[0.13,0.67,0.07,0.07],'FontSize',28,'Color','w','LineStyle','none','FontWeight','bold');
subplot(5,2,5); imagesc((abs(hdp3d(65:90,50:110,50,3)))); axis image; set(gca,'YTickLabels',[],'visible','off'); set(gca,'XTickLabels',[],'visible','off'); title('#3','FontSize',20); caxis([0 10]); cb = colorbar; cb.FontSize = 16; cb.Color = 'w';
annotation('textbox','String', 'R = 3','Position',[0.13,0.50,0.07,0.07],'FontSize',28,'Color','w','LineStyle','none','FontWeight','bold');
subplot(5,2,7); imagesc((abs(hdp3d(65:90,50:110,50,4)))); axis image; set(gca,'YTickLabels',[],'visible','off'); set(gca,'XTickLabels',[],'visible','off'); title('#4','FontSize',20); caxis([0 10]); cb = colorbar; cb.FontSize = 16; cb.Color = 'w';
annotation('textbox','String', 'R = 4','Position',[0.13,0.33,0.07,0.07],'FontSize',28,'Color','w','LineStyle','none','FontWeight','bold');
subplot(5,2,9); imagesc((abs(hdp3d(65:90,50:110,50,5)))); axis image; set(gca,'YTickLabels',[],'visible','off'); set(gca,'XTickLabels',[],'visible','off'); title('#5','FontSize',20); caxis([0 1]); cb = colorbar; cb.FontSize = 16; cb.Color = 'w';
annotation('textbox','String', 'R = 5','Position',[0.13,0.16,0.07,0.07],'FontSize',28,'Color','w','LineStyle','none','FontWeight','bold');
subplot(5,2,2); imagesc((abs(MRF_LRI5(65:90,50:110,1)))); axis image; set(gca,'YTickLabels',[],'visible','off'); set(gca,'XTickLabels',[],'visible','off'); title('\it{t_1}','FontSize',20); caxis([0 7]); cb = colorbar; cb.FontSize = 16; cb.Color = 'w';
annotation('textbox','String', '\it{t_0}','Position',[0.5926,0.84,0.07,0.07],'FontSize',28,'Color','w','LineStyle','none','FontWeight','bold');
subplot(5,2,4); imagesc((abs(MRF_LRI5(65:90,50:110,300)))); axis image; set(gca,'YTickLabels',[],'visible','off'); set(gca,'XTickLabels',[],'visible','off'); title('\it{t_1}','FontSize',20); caxis([0 7]); cb = colorbar; cb.FontSize = 16; cb.Color = 'w';
annotation('textbox','String', '\it{t_1}','Position',[0.5926,0.67,0.07,0.07],'FontSize',28,'Color','w','LineStyle','none','FontWeight','bold');
subplot(5,2,6); imagesc((abs(MRF_LRI5(65:90,50:110,600)))); axis image; set(gca,'YTickLabels',[],'visible','off'); set(gca,'XTickLabels',[],'visible','off'); title('\it{t_2}','FontSize',20); caxis([0 7]); cb = colorbar; cb.FontSize = 16; cb.Color = 'w';
annotation('textbox','String', '\it{t_2}','Position',[0.5926,0.50,0.07,0.07],'FontSize',28,'Color','w','LineStyle','none','FontWeight','bold');
subplot(5,2,8); imagesc((abs(MRF_LRI5(65:90,50:110,900)))); axis image; set(gca,'YTickLabels',[],'visible','off'); set(gca,'XTickLabels',[],'visible','off'); title('\it{t_3}','FontSize',20); caxis([0 7]); cb = colorbar; cb.FontSize = 16; cb.Color = 'w';
annotation('textbox','String', '\it{t_3}','Position',[0.5926,0.33,0.07,0.07],'FontSize',28,'Color','w','LineStyle','none','FontWeight','bold');
subplot(5,2,10); imagesc((abs(MRF_LRI5(65:90,50:110,1200)))); axis image; set(gca,'YTickLabels',[],'visible','off'); set(gca,'XTickLabels',[],'visible','off'); title('\it{t_0}','FontSize',20); caxis([0 7]); cb = colorbar; cb.FontSize = 16; cb.Color = 'w';
annotation('textbox','String', '\it{t_0}','Position',[0.5926,0.16,0.07,0.07],'FontSize',28,'Color','w','LineStyle','none','FontWeight','bold');
set(gcf,'color','k','InvertHardCopy','off'); colormap(gray());
annotation('textbox','String','(a)','FontSize',32,'LineStyle','none','Color','w','Position',[0.0123,0.8869,0.05,0.07])
annotation('textbox','String','(b)','FontSize',32,'LineStyle','none','Color','w','Position',[0.473,0.8869,0.05,0.07])

figure(3)
p1 = line([0 0],[0.4 1.8],'LineStyle','--','Color',[0.5 0.5 0.5],'LineWidth',3);
p2 = line([300 300],[0.4 1.8],'LineStyle','--','Color',[0.5 0.5 0.5],'LineWidth',3);
p3 = line([600 600],[0.4 1.6],'LineStyle','--','Color',[0.5 0.5 0.5],'LineWidth',3);
p4 = line([900 900],[0.4 1.6],'LineStyle','--','Color',[0.5 0.5 0.5],'LineWidth',3);
p5 = line([1200 1200],[0.4 1.8],'LineStyle','--','Color',[0.5 0.5 0.5],'LineWidth',3);
set(get(get(p1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(p2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(p3,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(p4,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(p5,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

text([0 300 600 900 1200]+10, 1.5*[1 1 1 1 1], {'\it{t_{0}}', '\it{t_{1}}', '\it{t_{2}}','\it{t_{3}}','\it{t_{0}}'},'FontSize',20,'Color',[0.5 0.5 0.5]);
shadedErrorBar([1:20:1200,1200],BSA_signal([1:20:1200,1200]),BSA_std([1:20:1200,1200]),'lineProps',{'color',cm(1,:),'LineWidth',2}); hold on;
shadedErrorBar([1:20:1200,1200],PL161_signal([1:20:1200,1200]),PL161_std([1:20:1200,1200]),'lineProps',{'color',cm(2,:),'LineWidth',2});
shadedErrorBar([1:20:1200,1200],MnCl2_signal([1:20:1200,1200]),MnCl2_std([1:20:1200,1200]),'lineProps',{'color',cm(3,:),'LineWidth',2});
plot([1:20:1200,1200],BSA_1DFT([1:20:1200,1200]),'o','LineWidth',2,'MarkerSize',8,'Color',cm(1,:));
plot([1:20:1200,1200],PL161_1DFT([1:20:1200,1200]),'o','LineWidth',2,'MarkerSize',8,'Color',cm(2,:));
plot([1:20:1200,1200],Water_1DFT([1:20:1200,1200]),'o','LineWidth',2,'MarkerSize',8,'Color',cm(3,:));
legend('BSA MT-MRF','PL161 MT-MRF','MnCl_{2} MT-MRF','BSA NPE','PL161 NPE','MnCl_{2} NPE','FontSize',16,'Orientation','horizontal'); legend boxoff; grid on; grid minor; ylim([0.4 1.8])
set(gca,'FontSize',20); xlabel('Pulse No.','FontSize',20); ylabel('Norm. Signal (a.u.)','FontSize',20)
set(gcf, 'Position',  [100, 100, 1920, 528])
annotation('textbox','String','(d)','FontSize',32,'LineStyle','none','Position',[0.0123,0.8869,0.05,0.07])

%% Reconstruction of in-vivo data for Figure 7.

close all; clear all; clc;

% This is the latest balanced data.
load('/home/dw16/Desktop/MRF_ihMT_Data/Mar20_Scans/hdp3d_ex.mat'); hdp3d = hdp3d_ex;

% Load dictionary.
load('/home/dw16/Desktop/MRF_ihMT_Data/Oblique_Phantom/Phantom_Dict_ex2.mat'); Ur = V(:,1:5);

ihMT = zeros(size(hdp3d,1),size(hdp3d,2),size(hdp3d,3)); MT = zeros(size(hdp3d,1),size(hdp3d,2),size(hdp3d,3));
for ii = 1:size(hdp3d,1)

    disp(['Slice number ', num2str(ii), ' out of ', num2str(size(hdp3d,1)), '.']);
    MRF_sing = squeeze(hdp3d(ii,:,:,1:5));
    mrf_aux = reshape(MRF_sing,[size(MRF_sing,1)*size(MRF_sing,2) size(MRF_sing,3)]); mrf_aux = mrf_aux';
    MRF_LRI5 = (Ur*mrf_aux)'; MRF_LRI5 = reshape(MRF_LRI5,[size(MRF_sing,1) size(MRF_sing,2) size(MRF_LRI5,2)]);
    
    % Take out global phase term.
    phase = repmat(angle(sum(MRF_LRI5,3)),[1 1 1200]);
    MRF_LRI5 = MRF_LRI5.*exp(-1i*phase);
    
    SB_av = (MRF_LRI5(:,:,600)+MRF_LRI5(:,:,1200))./2;
    TwoB_av = MRF_LRI5(:,:,300); ThreeB_av = MRF_LRI5(:,:,900);
   
    ihMT(ii,:,:) = real((TwoB_av-ThreeB_av)./(SB_av));
    MT(ii,:,:) = 1-(real((TwoB_av+ThreeB_av)./(2*SB_av)));

end

% Plot axial slice.
Slice = 150;

MRF_sing = squeeze(hdp3d(Slice,:,:,1:5));
mrf_aux = reshape(MRF_sing,[size(MRF_sing,1)*size(MRF_sing,2) size(MRF_sing,3)]); mrf_aux = mrf_aux';
MRF_LRI5 = (Ur*mrf_aux)'; MRF_LRI5 = reshape(MRF_LRI5,[size(MRF_sing,1) size(MRF_sing,2) size(MRF_LRI5,2)]);
brain_mask = squeeze(abs(MRF_LRI5(:,:,600))) > 1.5;

ihMT_map = squeeze(ihMT(Slice,:,:)); MT_map = squeeze(MT(Slice,:,:));
ihMT_map = ihMT_map.*brain_mask; MT_map = MT_map.*brain_mask;

figure(4)
subplot(2,3,5)
imagesc(imgaussfilt(ihMT_map(ceil(395/4):ceil(945/4),ceil(30/4):end-ceil(20/4))*100,0.7)); colormap(magma); caxis([0 5]); set(gca,'YTickLabels',[],'visible','off'); set(gca,'XTickLabels',[]); axis image 
cm = colorbar; cm.FontSize = 18; cm.Color = 'w'; title(cm,'ihMTR (%)','FontSize',18,'Color','w');
cm.Ticks = 0:1:4; cm.TickLabels = num2cell(0:1:4);
subplot(2,3,2)
imagesc(imgaussfilt(MT_map(ceil(395/4):ceil(945/4),ceil(30/4):end-ceil(20/4))*100,0.7)); colormap(magma); caxis([0 40]); set(gca,'YTickLabels',[],'visible','off'); set(gca,'XTickLabels',[]); axis image
cm = colorbar; cm.FontSize = 18; cm.Color = 'w'; title(cm,'MTR (%)','FontSize',18,'Color','w');
cm.Ticks = 0:10:30; cm.TickLabels = num2cell(0:10:30);

% Plot coronal slice.
Slice = 167;

MRF_sing = squeeze(hdp3d(:,Slice,:,1:5));
mrf_aux = reshape(MRF_sing,[size(MRF_sing,1)*size(MRF_sing,2) size(MRF_sing,3)]); mrf_aux = mrf_aux';
MRF_LRI5 = (Ur*mrf_aux)'; MRF_LRI5 = reshape(MRF_LRI5,[size(MRF_sing,1) size(MRF_sing,2) size(MRF_LRI5,2)]);
brain_mask = squeeze(abs(MRF_LRI5(:,:,600))) > 1.5;

ihMT_map = squeeze(ihMT(:,Slice,:)); MT_map = squeeze(MT(:,Slice,:));
ihMT_map = ihMT_map.*brain_mask; MT_map = MT_map.*brain_mask;

figure(4)
subplot(2,3,6)
imagesc(imgaussfilt(ihMT_map(ceil(350/4):ceil(785/4),ceil(25/4):ceil(465/4))*100,0.7)); colormap(magma); set(gca,'YTickLabels',[],'visible','off'); set(gca,'XTickLabels',[]); axis image 
caxis([0 5.5]); cm = colorbar; cm.FontSize = 18; cm.Color = 'w'; axis image; title(cm,'ihMTR (%)','FontSize',18,'Color','w');
subplot(2,3,3)
imagesc(imgaussfilt(MT_map(ceil(350/4):ceil(785/4),ceil(25/4):ceil(465/4))*100,0.7)); colormap(magma); set(gca,'YTickLabels',[],'visible','off'); set(gca,'XTickLabels',[]); axis image
caxis([0 42]); cm = colorbar; cm.FontSize = 18; cm.Color = 'w'; axis image; title(cm,'MTR (%)','FontSize',18,'Color','w');

% Plot sagittal slice.
Slice = 43;

MRF_sing = squeeze(hdp3d(:,:,Slice,1:5));
mrf_aux = reshape(MRF_sing,[size(MRF_sing,1)*size(MRF_sing,2) size(MRF_sing,3)]); mrf_aux = mrf_aux';
MRF_LRI5 = (Ur*mrf_aux)'; MRF_LRI5 = reshape(MRF_LRI5,[size(MRF_sing,1) size(MRF_sing,2) size(MRF_LRI5,2)]);
brain_mask = squeeze(abs(MRF_LRI5(:,:,600))) > 1.5;

ihMT_map = squeeze(ihMT(:,:,Slice)); MT_map = squeeze(MT(:,:,Slice));
ihMT_map = ihMT_map.*brain_mask; MT_map = MT_map.*brain_mask;

figure(4)
subplot(2,3,4)
imagesc(imgaussfilt(ihMT_map(ceil(355/4):ceil(810/4),ceil(400/4):ceil(930/4))*100,0.7)); caxis([0 5]); colormap(magma); set(gca,'YTickLabels',[],'visible','off'); set(gca,'XTickLabels',[]); axis image; 
cm = colorbar; cm.FontSize = 18; cm.Color = 'w'; title(cm,'ihMTR (%)','FontSize',18,'Color','w');
cm.Ticks = 0:1:4; cm.TickLabels = num2cell(0:1:4);
subplot(2,3,1)
imagesc(imgaussfilt(MT_map(ceil(355/4):ceil(810/4),ceil(400/4):ceil(930/4))*100,0.7)); caxis([0 40]); colormap(magma); set(gca,'YTickLabels',[],'visible','off'); set(gca,'XTickLabels',[]); axis image;
cm = colorbar; cm.FontSize = 18; cm.Color = 'w'; title(cm,'MTR (%)','FontSize',18,'Color','w');
cm.Ticks = 0:10:30; cm.TickLabels = num2cell(0:10:30);
set(gcf,'color','k','InvertHardCopy','off'); colormap(magma());

%% Plot singular images.

figure(5);
subplot(2,3,1)
imagesc(squeeze(abs(hdp3d(85:205,95:240,42,1)))); cb = colorbar; cb.FontSize = 16; cb.Color = 'w'; axis image
set(gca,'YTickLabels',[],'visible','off'); set(gca,'XTickLabels',[],'visible','off'); colormap(gray); caxis([0 200])
subplot(2,3,2)
imagesc(squeeze(abs(hdp3d(150,95:240,:,1)))); cb = colorbar; cb.FontSize = 16; cb.Color = 'w'; axis image
set(gca,'YTickLabels',[],'visible','off'); set(gca,'XTickLabels',[],'visible','off'); colormap(gray); caxis([0 200])
subplot(2,3,3)
imagesc(squeeze(abs(hdp3d(85:205,165,:,1)))); cb = colorbar; cb.FontSize = 16; cb.Color = 'w'; axis image
set(gca,'YTickLabels',[],'visible','off'); set(gca,'XTickLabels',[],'visible','off'); colormap(gray); caxis([0 200])
subplot(2,3,4)
imagesc(squeeze(abs(hdp3d(85:205,95:240,42,2)))); cb = colorbar; cb.Color = 'w'; cb.FontSize = 16; axis image
set(gca,'YTickLabels',[],'visible','off'); set(gca,'XTickLabels',[],'visible','off'); colormap(gray); caxis([0 20])
subplot(2,3,5)
imagesc(squeeze(abs(hdp3d(150,95:240,:,2)))); cb = colorbar; cb.FontSize = 16; cb.Color = 'w'; axis image
set(gca,'YTickLabels',[],'visible','off'); set(gca,'XTickLabels',[],'visible','off'); colormap(gray); caxis([0 20])
subplot(2,3,6)
imagesc(squeeze(abs(hdp3d(85:205,165,:,2)))); cb = colorbar; cb.FontSize = 16; cb.Color = 'w'; axis image
set(gca,'YTickLabels',[],'visible','off'); set(gca,'XTickLabels',[],'visible','off'); colormap(gray); caxis([0 20])
annotation('textbox','String', 'R = 1','Position',[0.03 0.82 0.1 0.1],'FontSize',36,'Color','w','LineStyle','none','FontWeight','bold');
annotation('textbox','String', 'R = 2','Position',[0.03 0.35 0.1 0.1],'FontSize',36,'Color','w','LineStyle','none','FontWeight','bold');
set(gcf,'color','k','InvertHardCopy','off'); colormap(gray());

figure(6);
subplot(2,3,1)
imagesc(squeeze(abs(hdp3d(85:205,95:240,42,3)))); cb = colorbar; cb.FontSize = 16; cb.Color = 'w'; axis image
set(gca,'YTickLabels',[],'visible','off'); set(gca,'XTickLabels',[],'visible','off'); colormap(gray); caxis([0 5])
subplot(2,3,2)
imagesc(squeeze(abs(hdp3d(150,95:240,:,3)))); cb = colorbar; cb.FontSize = 16; cb.Color = 'w'; axis image
set(gca,'YTickLabels',[],'visible','off'); set(gca,'XTickLabels',[],'visible','off'); colormap(gray); caxis([0 5])
subplot(2,3,3)
imagesc(squeeze(abs(hdp3d(85:205,165,:,3)))); cb = colorbar; cb.FontSize = 16; cb.Color = 'w'; axis image
set(gca,'YTickLabels',[],'visible','off'); set(gca,'XTickLabels',[],'visible','off'); colormap(gray); caxis([0 5])
subplot(2,3,4)
imagesc(squeeze(abs(hdp3d(85:205,95:240,42,4)))); cb = colorbar; cb.Color = 'w'; cb.FontSize = 16; axis image
set(gca,'YTickLabels',[],'visible','off'); set(gca,'XTickLabels',[],'visible','off'); colormap(gray); caxis([0 4])
subplot(2,3,5)
imagesc(squeeze(abs(hdp3d(150,95:240,:,4)))); cb = colorbar; cb.FontSize = 16; cb.Color = 'w'; axis image
set(gca,'YTickLabels',[],'visible','off'); set(gca,'XTickLabels',[],'visible','off'); colormap(gray); caxis([0 4])
subplot(2,3,6)
imagesc(squeeze(abs(hdp3d(85:205,165,:,4)))); cb = colorbar; cb.FontSize = 16; cb.Color = 'w'; axis image
set(gca,'YTickLabels',[],'visible','off'); set(gca,'XTickLabels',[],'visible','off'); colormap(gray); caxis([0 4])
annotation('textbox','String', 'R = 3','Position',[0.03 0.82 0.1 0.1],'FontSize',36,'Color','w','LineStyle','none','FontWeight','bold');
annotation('textbox','String', 'R = 4','Position',[0.03 0.35 0.1 0.1],'FontSize',36,'Color','w','LineStyle','none','FontWeight','bold');
set(gcf,'color','k','InvertHardCopy','off'); colormap(gray());

figure(7);
subplot(2,3,1)
imagesc(squeeze(abs(hdp3d(85:205,95:240,42,5)))); cb = colorbar; cb.FontSize = 16; cb.Color = 'w'; axis image
set(gca,'YTickLabels',[],'visible','off'); set(gca,'XTickLabels',[],'visible','off'); colormap(gray); caxis([0 2])
subplot(2,3,2)
imagesc(squeeze(abs(hdp3d(150,95:240,:,5)))); cb = colorbar; cb.FontSize = 16; cb.Color = 'w'; axis image
set(gca,'YTickLabels',[],'visible','off'); set(gca,'XTickLabels',[],'visible','off'); colormap(gray); caxis([0 2])
subplot(2,3,3)
imagesc(squeeze(abs(hdp3d(85:205,165,:,5)))); cb = colorbar; cb.FontSize = 16; cb.Color = 'w'; axis image
set(gca,'YTickLabels',[],'visible','off'); set(gca,'XTickLabels',[],'visible','off'); colormap(gray); caxis([0 2])
annotation('textbox','String', 'R = 5','Position',[0.03 0.82 0.1 0.1],'FontSize',36,'Color','w','LineStyle','none','FontWeight','bold');
set(gcf,'color','k','InvertHardCopy','off'); colormap(gray());
