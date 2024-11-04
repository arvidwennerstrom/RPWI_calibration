%% Choose directory and dates to include.
% Overall location of data
rootDir = 'C:\Users\arvwe\Onedrive - KTH\MEX\IRF'; % For my computer
% rootDir = 'C:\Users\arvidwen\Onedrive - KTH\MEX\IRF'; % For KTH computers
load([rootDir, '\Matlab\data\rpwi_data.mat'])


% Dates to include, on format 'YYMMDD'
date = '240819';
datasetDir = [rootDir '\datasets\20', date(1:2) '\', date(3:4), '\', date(5:6)];


% AC and DC differential mode: 0-3 (0 is default)
mode = 0;
[d_a,d_b,d_c] = mode_depending_LP_distance(mode,LP_delta_d);


data_handling_rules = struct('TM2voltage',true,'noise_filter',true,'density_mode_removal',true,'first_samples_removal',true,'voltage2Efield',true);
% data_handling_rules = struct('TM2voltage',false,'noise_filter',true,'density_mode_removal',true,'first_samples_removal',true,'voltage2Efield',false);


%% Read and handle data from .cdf's.
% Read and prepare HK-dat
HK10002_files = dir(fullfile(datasetDir, '\JUICE_LU_RPWI-PPTD-LWYHK10002*.cdf'));
HK10064_files = dir(fullfile(datasetDir, '\JUICE_LU_RPWI-PPTD-LWYHK10064*.cdf'));
HK10002_data = HK_data_prep(HK10002_files,10002);
HK10064_data = HK_data_prep(HK10064_files,10064);

% Read and prepare LP-dat
LP_files = dir(fullfile(datasetDir, '\JUICE_L1a_RPWI-LP-SID1_RICH_DE763_SNAP*.cdf'));
for i=1:length(LP_files)
    fullPath = fullfile(LP_files(i).folder, LP_files(i).name);
    LP_file = dataobj(fullPath);
    Epoch_data = LP_file.data.Epoch.data();
    convertedEpoch = EpochTT(Epoch_data);

    % Reading raw data
    LP_TM = LP_file;
    adcs = double(LP_file.data.LP_DATA.data);
    
    
    if data_handling_rules.TM2voltage
    % Calibration using Lars' coefficients.
        adcs = TM_2_voltage(adcs);
        if data_handling_rules.voltage2Efield
            % Convert voltages to Electric field strength between probes.         
            adcs = voltage_2_Efield(adcs,[d_a(4) d_b(4) d_c(4)]);
        end
    end


    
    
         
    if data_handling_rules.density_mode_removal
        % Remove data where LPs are in density mode @ 20240126.
        [adcs,convertedEpoch] = density_mode_removal(adcs,convertedEpoch,LP_files(i).name);
    end


    if data_handling_rules.first_samples_removal
        % Remove the first samples of the measurement.  
        cut_width = 32; % # of samples to remove whenever needed.
        if i == 1
            adcs = adcs(cut_width+1:end,:);
            convertedEpoch = convertedEpoch(cut_width+1:end);
        end
    end

    
    % Combine the data
    m=length(convertedEpoch);
    data_differentials = irf.ts_scalar(convertedEpoch, adcs(:,1:3));
    data_single = irf.ts_scalar(convertedEpoch, adcs(:,4));
    data_all = irf.ts_scalar(Epoch_data, LP_TM.data.LP_DATA.data);

    if i == 1
        LP_data_differentials = data_differentials;
        LP_data_single = data_single;
        LP_data_All = data_all;
    else
        LP_data_differentials = combine(LP_data_differentials, data_differentials);
        LP_data_single = combine(LP_data_single, data_single);
        LP_data_All = combine(LP_data_All,data_all);
    end
end


%% Maybe add the v7.3
save([rootDir, '\Matlab\data\data_', date],'HK10002_files','HK10002_data','HK10064_files','HK10064_data','LP_files','LP_data_All','LP_data_differentials','LP_data_single','LP_TM','data_handling_rules')


%% Functions
function [d_a,d_b,d_c] = mode_depending_LP_distance(mode,LP_distances)
    if mode == 0
        d_a = LP_distances.d_12; d_b = LP_distances.d_23; d_c = LP_distances.d_34;
    elseif mode == 1
        d_a = LP_distances.d_14; d_b = LP_distances.d_13; d_c = LP_distances.d_34;
    elseif mode == 2
        d_a = LP_distances.d_12; d_b = LP_distances.d_24; d_c = LP_distances.d_14;
    elseif mode == 3
        d_a = LP_distances.d_14; d_b = LP_distances.d_14; d_c = LP_distances.d_14;
    end
end

function [data] = TM_2_voltage(data)
    % Coefficients for channel # in order: [A B C D] 
    calibration_coefficients = [5.15*1e-6, 4.97*1e-6, 5.07*1e-6, 9.94*1e-5];
    for channel = 1:width(data)
        data(:,channel) = calibration_coefficients(channel)*(data(:,channel));
    end
end



function [data,time] = density_mode_removal(data,time,filename)
    % Seq. 14: Index 1379810:1927040 in *T101713, 
    % Seq. 17: Index 2371426:end in *T101713 & 1:399969 in *T112039,
    % Seq. 18 and 20: Index 1737922:end in *T112039 
    if filename == 'JUICE_L1a_RPWI-LP-SID1_RICH_DE763_SNAP_20240126T101713_V03.cdf'
        data = [data(1:1379809,:); data(1927041:2371425,:)];
        time = [time(1:1379809,:); time(1927041:2371425,:)];    
    elseif filename == 'JUICE_L1a_RPWI-LP-SID1_RICH_DE763_SNAP_20240126T112039_V02.cdf'
        data = data(399970:1737921,:);
        time = time(399970:1737921,:);
    end
end

function [data] = voltage_2_Efield(data,distance)
        % Get E-field strength using data (voltage) and LP distances.
        % adcs, since calibration above, contain values in V, distances 
        % "d_a/b/c" in meters. Adding "1e3" makes new unit of adcs
        % mV/m.
        data(:,1) = 1e3*data(:,1)/distance(1);
        data(:,2) = 1e3*data(:,2)/distance(2);
        data(:,3) = 1e3*data(:,3)/distance(3);


end

