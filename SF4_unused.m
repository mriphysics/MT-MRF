%%% ts-ihMT and ss-ihMT signal simulation and comparison.

close all; clear all; clc;

n_increments = 50;
T1f_range = linspace(0.2,4,n_increments);
T1s_range = linspace(0.1,1,n_increments);
T1d_range = linspace(2e-3,25e-3,n_increments);
T2f_range = linspace(50e-3,500e-3,n_increments);
T2s_range = linspace(5e-6,20e-6,n_increments);
delta_range = linspace(0,1,n_increments);
k_range = linspace(40,100,n_increments);
M0s_range = linspace(0,0.25,n_increments);

param_ranges = [T1f_range ; T1s_range ; T1d_range ; T2f_range ; T2s_range ; delta_range ; k_range ; M0s_range];

ihMTR_ts = zeros(n_increments,8); B1max = zeros(n_increments,8);

for pp = 1:8
    
    var_param = param_ranges(pp,:);

    for ii = 1:n_increments
        
        if pp == 1
            T12x = [var_param(ii) 1000*1e-3 6.5e-3 80e-3 12.5e-6];
            f = 0.65; K = 65; M0b_Varma = (7.3*(1/T12x(1)))/K; M0b = M0b_Varma/(1 + M0b_Varma);
            M0f = 1-M0b; M0 = [M0f M0b*(1-f) M0b*f];
        elseif pp == 2
            T12x = [650e-3 var_param(ii) 6.5e-3 80e-3 12.5e-6];
            f = 0.65; K = 65; M0b_Varma = (7.3*(1/T12x(1)))/K; M0b = M0b_Varma/(1 + M0b_Varma);
            M0f = 1-M0b; M0 = [M0f M0b*(1-f) M0b*f];
        elseif pp == 3
            T12x = [650e-3 1000*1e-3 var_param(ii) 80e-3 12.5e-6];
            f = 0.65; K = 65; M0b_Varma = (7.3*(1/T12x(1)))/K; M0b = M0b_Varma/(1 + M0b_Varma);
            M0f = 1-M0b; M0 = [M0f M0b*(1-f) M0b*f];
        elseif pp == 4
            T12x = [650e-3 1000*1e-3 6.5e-3 var_param(ii) 12.5e-6];
            f = 0.65; K = 65; M0b_Varma = (7.3*(1/T12x(1)))/K; M0b = M0b_Varma/(1 + M0b_Varma);
            M0f = 1-M0b; M0 = [M0f M0b*(1-f) M0b*f];
        elseif pp == 5
            T12x = [650e-3 1000*1e-3 6.5e-3 80e-3 var_param(ii)];
            f = 0.65; K = 65; M0b_Varma = (7.3*(1/T12x(1)))/K; M0b = M0b_Varma/(1 + M0b_Varma);
            M0f = 1-M0b; M0 = [M0f M0b*(1-f) M0b*f];
        elseif pp == 6
            T12x = [650e-3 1000*1e-3 6.5e-3 80e-3 12.5e-6];
            f = var_param(ii); K = 65; M0b_Varma = (7.3*(1/T12x(1)))/K; M0b = M0b_Varma/(1 + M0b_Varma);
            M0f = 1-M0b; M0 = [M0f M0b*(1-f) M0b*f];
        elseif pp == 7
            T12x = [650e-3 1000*1e-3 6.5e-3 80e-3 12.5e-6];
            f = 0.65; K = var_param(ii); M0b_Varma = (7.3*(1/T12x(1)))/K; M0b = M0b_Varma/(1 + M0b_Varma);
            M0f = 1-M0b; M0 = [M0f M0b*(1-f) M0b*f];
        elseif pp == 8
            T12x = [650e-3 1000*1e-3 6.5e-3 80e-3 12.5e-6];
            f = 0.65; K = 65; M0b = var_param(ii);
            M0f = 1-M0b; M0 = [M0f M0b*(1-f) M0b*f];
        end
        
        flips = d2r(29.51); TR = 5.33e-3; Dur = 2.51e-3; Delta = 8058.48; n1B = 300; nMB = 300;
        TBW = 2; B1rms = 4; dphi = 0; Var_Params = ([M0, T12x, K]);
        
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
        
        ihMTR_ts(pp,ii) = abs((Mss_ts(nMB+n1B+nMB)-Mss_ts(nMB))/((Mss_ts(nMB+n1B)+Mss_ts(nMB+n1B+nMB+n1B))/2)*100);
        
    end 

end

cm = lines(7);
figure(1)
subplot(2,4,1);
plot(linspace(0.2,4,n_increments),ihMTR_ts(1,:),'-','Color',cm(1,:),'LineWidth',3); xlabel('T_{1}^{f} (s)'), ylabel('ihMTR (%)'); set(gca,'FontSize',18); grid on; grid minor; xlim([0.2 4]); ylim([0 15]);
subplot(2,4,2);
plot(linspace(0.1,1,n_increments),ihMTR_ts(2,:),'-','Color',cm(2,:),'LineWidth',3); xlabel('T_{1Z}^{s} (s)'), ylabel('ihMTR (%)'); set(gca,'FontSize',18); grid on; grid minor; xlim([0.1 1]); ylim([0 15]);
subplot(2,4,3);
plot(linspace(2,25,n_increments),ihMTR_ts(3,:),'-','Color',cm(3,:),'LineWidth',3); xlabel('T_{1D}^{s} (ms)'), ylabel('ihMTR (%)'); set(gca,'FontSize',18); grid on; grid minor; xlim([2 25]); ylim([0 15]);
subplot(2,4,4);
plot(linspace(50,500,n_increments),ihMTR_ts(4,:),'-','Color',cm(4,:),'LineWidth',3); xlabel('T_{2}^{f} (ms)'), ylabel('ihMTR (%)'); set(gca,'FontSize',18); grid on; grid minor; xlim([50 500]); ylim([0 15]);
subplot(2,4,5);
plot(linspace(5,20,n_increments),ihMTR_ts(5,:),'-','Color',cm(5,:),'LineWidth',3); xlabel('T_{2}^{s} (\mus)'), ylabel('ihMTR (%)'); set(gca,'FontSize',18); grid on; grid minor; xlim([5 20]); ylim([0 15]);
subplot(2,4,6);
plot(linspace(0,1,n_increments),ihMTR_ts(6,:),'-','Color',cm(6,:),'LineWidth',3); xlabel('\delta'), ylabel('ihMTR (%)'); set(gca,'FontSize',18); grid on; grid minor; xlim([0 1]); ylim([0 15]);
subplot(2,4,7);
plot(linspace(40,100,n_increments),ihMTR_ts(7,:),'-','Color',cm(7,:),'LineWidth',3); xlabel('\it{k} (s^{-1})'), ylabel('ihMTR (%)'); set(gca,'FontSize',18); grid on; grid minor; xlim([40 100]); ylim([0 15]);
subplot(2,4,8);
plot(linspace(0,0.25,n_increments),ihMTR_ts(8,:),'-','Color','k','LineWidth',3); xlabel('M_{0}^{s}'), ylabel('ihMTR (%)'); set(gca,'FontSize',18); grid on; grid minor; xlim([0 0.25]); ylim([0 15]);
