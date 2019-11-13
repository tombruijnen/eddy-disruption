%% Make plots for the girf computation
figure('units','normalized','outerposition',[0 0 1 1]);
subplot(221);plot(f_storage,abs(girf_zeroth(:,:)),'LineWidth',2);
ylabel('Magnitude [a.u]');xlabel('Frequency [Hz]');title('|GIRF_{0}|');
legend('X','Y','Z','orientation','horizontal','location','north')
grid on;set(gca,'FontSize',14,'LineWidth',2);axis([])

subplot(222);plot(f_storage,angle(girf_zeroth(:,:)),'LineWidth',2);
ylabel('Argument [rad]');xlabel('Frequency [Hz]');title('Arg(GIRF_{0})');
legend('X','Y','Z','orientation','horizontal','location','north')
grid on;set(gca,'FontSize',14,'LineWidth',2);axis([])

subplot(223);plot(f_storage,abs(girf_first),'LineWidth',2);
ylabel('Magnitude [a.u.]');xlabel('Frequency [Hz]');title('|GIRF_1|')
legend('X','Y','Z','orientation','horizontal','location','north')
grid on;set(gca,'FontSize',14,'LineWidth',2);

subplot(224);plot(imag(girf_first),-real(girf_first),'LineWidth',2);
ylabel('Real [a.u.]');xlabel('Imag [a.u.]');title('GIRF_1')
legend('X','Y','Z','orientation','horizontal','location','south')
grid on;set(gca,'FontSize',14,'LineWidth',2);

set(gcf,'Color','w')