clear all
%% Choose directory and dates to include.
% Overall location of data
rootDir = 'C:\Users\arvwe\Onedrive - KTH\MEX\IRF'; % For my computer
% rootDir = 'C:\Users\arvidwen\Onedrive - KTH\MEX\IRF'; % For KTH computers

date = '240820';

load([rootDir, '\Matlab\data\data_', date])
%% Plotting.
figure('Name',['LP data: ' date],'NumberTitle','off')
subplot(2,1,1); irf_plot(LP_data_differentials,'.')
legend('P12', 'P23', 'P34')
if data_handling_rules.TM2voltage
    if data_handling_rules.voltage2Efield
        ylabel('mV/m')
    else
        ylabel('V')
    end
else
    ylabel('TM')
end
xlabel('Time')
title('Probe differentials');
subplot(2,1,2); irf_plot(LP_data_single, '.')
legend('P4')
if data_handling_rules.TM2voltage
    ylabel('V')
else
    ylabel('TM')
end
xlabel('Time')
title('Single ended probe (P4) voltage');

figure('Name',['HK data: ' date],'NumberTitle','off')
subplot(2,1,1); irf_plot(HK10002_data, '*')
legend('P1', 'P2', 'P3', 'P4')
ylabel('°C')
xlabel('Time')
title('LWYHK10002: Temperatures');

subplot(2,1,2); irf_plot(HK10064_data, '*')
legend('LP')
ylabel('°C')
xlabel('Time')
title('LWYHK10064: Temperatures');

E_field_aligned = align_Efield(rootDir,LP_data_differentials.data,0);
figure()
plot(1:length(E_field_aligned),E_field_aligned,'.')
legend('x-axis','y-axis','z-axis')
ylabel('mV/m')
title("Electric field in SC primary axes")