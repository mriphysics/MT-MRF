%% Script to simulate transient and steady-state signals for Figure 2. Daniel West 2020.

close all; clear all; clc;

T12x = [650*1e-3 1000*1e-3 6.5e-3 80e-3 12.5e-6];
delta = 0.65; k = 65; M0s_Varma = (7.3*(1/T12x(1)))/k; M0s = M0s_Varma/(1 + M0s_Varma);
M0f = 1-M0s; M0 = [M0f M0s*(1-delta) M0s*delta];

%% Experimental acquisition scheme.

flips = d2r(29.51); TR = 5.33e-3; Dur = 2.51e-3; Delta = 8058.48; n1B = 300; nMB = 300;
TBW = 2; B1rms = 4; dphi = 0; Var_Params = ([M0, T12x, k]);

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

% Pre-compute lineshape for different T2b - removed loop because T2b fixed.
df_1B = [0 0 0]; df_2B = [0 0 Delta]; df_3B = [-Delta 0 Delta];
[G_3B,wloc_3B] = SuperLorentzian_LSint(T12x(5),df_3B);
[G_2B,wloc_2B] = SuperLorentzian_LSint(T12x(5),df_2B);
[G_1B,wloc_1B] = SuperLorentzian_LSint(T12x(5),df_1B);

Mss_ts = Dictionary_function_CSS(flips,TR,Dur,dphi,Var_Params(4:6),Var_Params(7:8),Var_Params(1:3),Var_Params(9),Delta,TBW,nMB,n1B,b1sqrd,G_1B,G_2B,G_3B,wloc_1B,wloc_2B,wloc_3B);
ihMTR_ts = abs((Mss_ts(nMB+n1B+nMB)-Mss_ts(nMB))/((Mss_ts(nMB+n1B)+Mss_ts(nMB+n1B+nMB+n1B))/2)*100);

%% Optimal acquisition scheme.

flips = d2r(29.51); TR = 5.33e-3; Dur = 2.51e-3; Delta = 8058.48; n1B = 152; nMB = 141;
TBW = 2; B1rms = 4; dphi = 0; Var_Params = ([M0, T12x, k]);

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

% Pre-compute lineshape for different T2b - removed loop because T2b fixed.
df_1B = [0 0 0]; df_2B = [0 0 Delta]; df_3B = [-Delta 0 Delta];
[G_3B,wloc_3B] = SuperLorentzian_LSint(T12x(5),df_3B);
[G_2B,wloc_2B] = SuperLorentzian_LSint(T12x(5),df_2B);
[G_1B,wloc_1B] = SuperLorentzian_LSint(T12x(5),df_1B);

Mss_op = Dictionary_function_CSS(flips,TR,Dur,dphi,Var_Params(4:6),Var_Params(7:8),Var_Params(1:3),Var_Params(9),Delta,TBW,nMB,n1B,b1sqrd,G_1B,G_2B,G_3B,wloc_1B,wloc_2B,wloc_3B);
ihMTR_op = abs((Mss_op(nMB+n1B+nMB)-Mss_op(nMB))/((Mss_op(nMB+n1B)+Mss_op(nMB+n1B+nMB+n1B))/2)*100);

%% Forbidden acquisition scheme.

flips = d2r(29.51); TR = 5.33e-3; Dur = 2.51e-3; Delta = 8058.48; n1B = 350; nMB = 100;
TBW = 2; B1rms = 4; dphi = 0; Var_Params = ([M0, T12x, k]);

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

% Pre-compute lineshape for different T2b - removed loop because T2b fixed.
df_1B = [0 0 0]; df_2B = [0 0 Delta]; df_3B = [-Delta 0 Delta];
[G_3B,wloc_3B] = SuperLorentzian_LSint(T12x(5),df_3B);
[G_2B,wloc_2B] = SuperLorentzian_LSint(T12x(5),df_2B);
[G_1B,wloc_1B] = SuperLorentzian_LSint(T12x(5),df_1B);

Mss_fo = Dictionary_function_CSS(flips,TR,Dur,dphi,Var_Params(4:6),Var_Params(7:8),Var_Params(1:3),Var_Params(9),Delta,TBW,nMB,n1B,b1sqrd,G_1B,G_2B,G_3B,wloc_1B,wloc_2B,wloc_3B);
ihMTR_fo = abs((Mss_fo(nMB+n1B+nMB)-Mss_fo(nMB))/((Mss_fo(nMB+n1B)+Mss_fo(nMB+n1B+nMB+n1B))/2)*100);

%% ss-ihMT simulation

flips = d2r(30); TR = 5e-3; Dur = 2.3e-3; Delta = 8000; n1B = 0; nMB = 1200; % Trick to get the sequence to simulate the full 2B and 3B periods.
TBW = 2; B1rms = 4; dphi = 0; Var_Params = ([M0, T12x, k]);

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
[G_3B,wloc_3B] = SuperLorentzian_LSint(T12x(5),df_3B);
[G_2B,wloc_2B] = SuperLorentzian_LSint(T12x(5),df_2B);
[G_1B,wloc_1B] = SuperLorentzian_LSint(T12x(5),df_1B);

% This yields the two steady-state magnetisation using each pulse type.
Mss_Total_ss = Dictionary_function_CSS(flips,TR,Dur,dphi,Var_Params(4:6),Var_Params(7:8),Var_Params(1:3),Var_Params(9),Delta,TBW,nMB,n1B,b1sqrd,G_1B,G_2B,G_3B,wloc_1B,wloc_2B,wloc_3B);
Mss_SB = Dictionary_function_CSS(flips,TR,Dur,dphi,Var_Params(4:6),Var_Params(7:8),Var_Params(1:3),Var_Params(9),Delta,TBW,0,1200,b1sqrd,G_1B,G_2B,G_3B,wloc_1B,wloc_2B,wloc_3B);

% Extract the final values corresponding to a sequence of 1200 2B and 1200 3B pulses.
Mss_2Bss = Mss_Total_ss(1200)*ones(1,1200);
Mss_3Bss = Mss_Total_ss(2400)*ones(1,1200);

ihMTR_ss = (abs(Mss_2Bss(end)-Mss_3Bss(end))/abs(Mss_SB(end)))*100;

%% Plot results.

Line_2B = abs(Mss_ts(300)*ones(1,1200));
Line_3B = abs(Mss_ts(900)*ones(1,1200));

Dihmt_ss = abs(Mss_2Bss(end))-abs(Mss_3Bss(end));
Dihmt_ts = abs(Mss_ts(300))-abs(Mss_ts(900));

figure(1);
cm = colormap(lines(5));
plot(abs(Mss_op),'LineStyle','-','LineWidth',3,'Color',[0.4940    0.1840    0.5560]); hold on;
plot(abs(Mss_fo),'LineStyle','-','LineWidth',3,'Color',[0.4660    0.6740    0.1880]);
plot(abs(Mss_ts),'LineStyle','-','LineWidth',3,'Color',[0.9290    0.6940    0.1250]); hold on;
plot(abs(Mss_3Bss),'LineStyle','-','LineWidth',3,'Color',cm(2,:));
plot(abs(Mss_2Bss),'LineStyle','-','LineWidth',3,'Color',cm(1,:));
plot(Line_2B,'LineStyle','--','LineWidth',2,'Color',[0.5 0.5 0.5]);
plot(Line_3B,'LineStyle','--','LineWidth',2,'Color',[0.5 0.5 0.5]);
xlabel('Pulse No.','FontSize',20); ylabel('Signal (Fraction of M_{0})','FontSize',20);
get(gca,'XTick'); set(gca, 'FontSize',18); get(gca,'YTick'); set(gca, 'FontSize',18); grid on; grid minor
ll = legend('Maximum Contrast','Forbidden','Scan','ss-ihMT 3B','ss-ihMT 2B'); ll.FontSize = 20; legend boxoff; ylim([0.06 0.16]); xlim([1 1200]); ll.Orientation = 'horizontal';

annotation('doublearrow',[0.8 0.8],[0.433 0.473],'HeadStyle','plain','LineWidth',1)
annotation('doublearrow',[0.8 0.8],[0.27 0.348],'HeadStyle','plain','LineWidth',1)

text(1050, 0.102,'\DeltaihMT = 0.006','FontSize',15);
text(1050, 0.084,'\DeltaihMT = 0.011','FontSize',15);
