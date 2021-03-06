%% Make plots for the signal simulations
% Signal intensity off-resonance landscapes
figure('units','normalized','outerposition',[0 0 1 1])
subplot(221);imshow(abs(signal),[]);colormap hot;colorbar;title('No eddy current');xlabel(['\DeltaB_{0} = [',num2str(sim.freq(1)),'; ',num2str(sim.freq(end)),'] Hz']);ylabel(['# RF pulse = [1; ',num2str(sim.N),']']);set(gca,'Fontsize',16);
subplot(222);imshow(abs(signal_lin),[]);colormap hot;colorbar;title('Eddy + LIN');xlabel(['\DeltaB_{0} = [',num2str(sim.freq(1)),'; ',num2str(sim.freq(end)),'] Hz']);ylabel(['# RF pulse = [1; ',num2str(sim.N),']']);set(gca,'Fontsize',16)
subplot(223);imshow(abs(signal_rnd),[]);colormap hot;colorbar;title('Eddy + RND');xlabel(['\DeltaB_{0} = [',num2str(sim.freq(1)),'; ',num2str(sim.freq(end)),'] Hz']);ylabel(['# RF pulse = [1; ',num2str(sim.N),']']);set(gca,'Fontsize',16)
subplot(224);imshow(abs(signal_rad),[]);colormap hot;colorbar;title('Eddy + GA-RAD');xlabel(['\DeltaB_{0} = [',num2str(sim.freq(1)),'; ',num2str(sim.freq(end)),'] Hz']);ylabel(['# RF pulse = [1; ',num2str(sim.N),']']);set(gca,'Fontsize',16)

% Signal intensity profiles with standard deviation
incr = 1;
rang = 100;
axs = [sim.freq(1) sim.freq(end) 0 1.5*max(matrix_to_vec(abs(signal(end-rang:end,:))))];
figure('units','normalized','outerposition',[0 0 1 1])
subplot(221);
shadedErrorBar(sim.freq(1:incr:end),abs(mean(signal(end-rang:end,1:incr:end),1)),...
    abs(std(signal(end-rang:end,1:incr:end),0,1)),'lineprops',{'Color','b','LineWidth',3});hold on
ylabel('|Signal Intensity|');grid on;xlabel('\DeltaB_{0} [Hz]');set(gca,'FontSize',16,'LineWidth',3);axis(axs);title('No eddy current')

subplot(222);
shadedErrorBar(sim.freq(1:incr:end),abs(mean(signal_lin(end-rang:end,1:incr:end),1)),...
    abs(std(signal_lin(end-rang:end,1:incr:end),0,1)),'lineprops',{'Color','b','LineWidth',3});hold on
ylabel('|Signal Intensity|');grid on;xlabel('\DeltaB_{0} [Hz]');set(gca,'FontSize',16,'LineWidth',3);axis(axs);title('Eddy + LIN')

subplot(223);
shadedErrorBar(sim.freq(1:incr:end),abs(mean(signal_rnd(end-rang:end,1:incr:end),1)),...
    abs(std(signal_rnd(end-rang:end,1:incr:end),0,1)),'lineprops',{'Color','b','LineWidth',3});hold on
ylabel('|Signal Intensity|');grid on;xlabel('\DeltaB_{0} [Hz]');set(gca,'FontSize',16,'LineWidth',3);axis(axs);title('Eddy + RND')

subplot(224);
shadedErrorBar(sim.freq(1:incr:end),abs(mean(signal_rad(end-rang:end,1:incr:end),1)),...
    abs(std(signal_rad(end-rang:end,1:incr:end),0,1)),'lineprops',{'Color','b','LineWidth',3});hold on
ylabel('|Signal Intensity|');grid on;xlabel('\DeltaB_{0} [Hz]');set(gca,'FontSize',16,'LineWidth',3);axis(axs);title('Eddy + GA-RAD')
