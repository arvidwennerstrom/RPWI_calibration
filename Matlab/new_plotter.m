clear all
%% Choose directory and dates to include.
% Overall location of data
rootDir = 'C:\Users\arvwe\Onedrive - KTH\MEX\IRF'; % For my computer
% rootDir = 'C:\Users\arvidwen\Onedrive - KTH\MEX\IRF'; % For KTH computers


% Dates to include, on format 'YYMMDD'
date = '240820';


load([rootDir, '\P04_20240820T200000.mat']);
lpstruc1 = load([rootDir, '\sweep_data\lpstruc_20240820_1_20241002.mat']);
lpstruc2 = load([rootDir, '\sweep_data\lpstruc_20240820_2_20241002.mat']);
% load(rootDir, '\sweep_data\lpstruc_20240820_2_20241002.mat')
%% Plotting
t = lpstruc2.lpstruc.t_sweep;
Vi = lpstruc2.lpstruc.Vi;
Ni = lpstruc2.lpstruc.Ni;
Usc = lpstruc2.lpstruc.Usc

figure()
plot(t,-Usc);
% hold on
% plot(t,Ni)