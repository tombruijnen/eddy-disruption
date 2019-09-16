%% Make plots for the brodsky and duyn measurements
figure('units','normalized','outerposition',[0 0 1 1]);
subplot(221);plot(10^3*labels.t,10^3 * labels.G(:,1),'LineWidth',2);
ylabel('G [mT/m]');xlabel('Time [ms]');title('Input waveform');
grid on;set(gca,'FontSize',14,'LineWidth',2);

subplot(222);plot(10^3*labels.t_adc,squeeze(kdata(:,1,:,1)),'LineWidth',2);
ylabel('\Delta\phi [rad]');xlabel('Time [ms]');title('Raw kdata')
grid on;set(gca,'FontSize',14,'LineWidth',2);

subplot(223);plot(10^3*labels.t_adc,squeeze(zeroth_order(:,1,1)),'LineWidth',2);
ylabel('\DeltaB_0 [Hz]');xlabel('Time [ms]');title('Brodsky analysis')
grid on;set(gca,'FontSize',14,'LineWidth',2)

subplot(224);plot(10^3*labels.t_adc,10^3*squeeze(first_order(:,1,1)),'LineWidth',2);
ylabel('G [mT/m] [rad]');xlabel('Time [ms]');title('Duyn analysis')
grid on;set(gca,'FontSize',14,'LineWidth',2)