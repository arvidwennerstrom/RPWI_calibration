function [data_all] = HK_data_prep(my_files,HK_type)
 
    data_all = [];
    
    for i=1:length(my_files)
        fullPath = fullfile(my_files(i).folder, my_files(i).name);
        HK_file = dataobj(fullPath);
        Epoch_data = HK_file.data.Epoch.data();
        convertedEpoch = EpochTT(Epoch_data);
        
        % Create different data based on type of HK file. 
        if HK_type == 10002
            adc_a = HK_file.data.LWT0449D_CALIBRATED.data;
            adc_b = HK_file.data.LWT0449E_CALIBRATED.data;
            adc_c = HK_file.data.LWT0449F_CALIBRATED.data;
            adc_d = HK_file.data.LWT044A0_CALIBRATED.data;
            
            n=1;
            m=length(HK_file.data.Epoch.data());
            data = irf.ts_scalar(convertedEpoch(n:m), [adc_a(n:m), adc_b(n:m), adc_c(n:m), adc_d(n:m)]);
        
        elseif HK_type == 10064
            adc_a = HK_file.data.LWT04567.data;

            n=1;
            m=length(HK_file.data.Epoch.data());
            data = irf.ts_scalar(convertedEpoch(n:m), [adc_a(n:m)]);
        end

        if i == 1
            data_all = data;
        else
            data_all = combine(data_all, data);
        end
    end
end