figure(3), clf
randFID = round(rand(1,20)*1000)

subplot(3,1,1)
plot(abs(dissolved_data(1:32,randFID)))
ylabel('Magnitude')
title('20 Random Dissolved FIDs')

hold on
yl = ylim
% plot([64 64],[yl(1) yl(2)],'--k')

subplot(3,1,2)
plot(real(dissolved_data(1:32,randFID)))
ylabel('Real')

hold on
yl = ylim
% plot([64 64],[yl(1) yl(2)],'--k')

subplot(3,1,3)
plot(imag(dissolved_data(1:32,randFID)))
ylabel('Imaginary')
xlabel('Points in FID')

hold on
yl = ylim
% plot([64 64],[yl(1) yl(2)],'--k')

%%
figure(2), clf
randFID = round(rand(1,20)*1000)

subplot(3,1,1)
plot(abs(gas_data(:,randFID))./max(abs(gas_data(:,randFID))))
ylabel('Magnitude')
title('20 Random Gas FIDs')

hold on
yl = ylim
set(gca, 'YScale', 'log')
% plot([64 64],[yl(1) yl(2)],'--k')

subplot(3,1,2)
plot(real(gas_data(:,randFID)))
ylabel('Real')
set(gca, 'YScale', 'log')

hold on
yl = ylim
% plot([64 64],[yl(1) yl(2)],'--k')

subplot(3,1,3)
plot(imag(gas_data(:,randFID)))
ylabel('Imaginary')
xlabel('Points in FID')
set(gca, 'YScale', 'log')
hold on
yl = ylim
% plot([64 64],[yl(1) yl(2)],'--k')
%%
load('D:\OneDrive\Documents\Duke\CIVM Research\keyhole\GE_002065.mat')
GE_gas_data = gas_data;
[dissolved_data, gas_data,  dwell_time] = prep_data_for_keyhole(all_keyhole{1,22});

figure(2), clf
randFID = round(rand(1,20)*1000)

subplot(2,1,1)
plot(abs(GE_gas_data(:,randFID))./max(abs(GE_gas_data(:,randFID))))
ylabel('Magnitude')
title('Normalized GE Gas FIDs')
% set(gca, 'YScale', 'log')

subplot(2,1,2)
plot(abs(gas_data(1:64,randFID))./max(abs(gas_data(:,randFID))))
ylabel('Magnitude')
title('Normalized Siemens Gas FIDs')
% set(gca, 'YScale', 'log')