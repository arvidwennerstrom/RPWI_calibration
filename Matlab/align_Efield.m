function [E_field_aligned] = align_Efield(rootDir,E_field,mode)
%     E_field is electric field strength for differentials a,b and c [mV/m]

    load([rootDir, '\Matlab\data\rpwi_data.mat'])
    if mode == 0
        d_a = LP_delta_d.d_12; d_b = LP_delta_d.d_23; d_c = LP_delta_d.d_34;
    elseif mode == 1
        d_a = LP_delta_d.d_14; d_b = LP_delta_d.d_13; d_c = LP_delta_d.d_34;
    elseif mode == 2
        d_a = LP_delta_d.d_12; d_b = LP_delta_d.d_24; d_c = LP_delta_d.d_14;
    elseif mode == 3
        d_a = LP_delta_d.d_14; d_b = LP_delta_d.d_14; d_c = LP_delta_d.d_14;
    end
    
%   M is confirmed by "The Radio & Plasma Wave Investigation (RPWI) 
%   for the JUpiter ICy moons Explorer (JUICE)", page 34
%   M^-1 = [0.1852 0.1923 0.1917; 0.1320 -0.0112 -0.0322; 0.0 0.1398 0.0];
    M = [d_a(1:3); 
        d_b(1:3);
        d_c(1:3)]; % [m]
    
    % d_a(4)
%   Converted from [mV/m] to [mV] by multiplying with "d_a/b/c".
    deltaU = [E_field(:,1)*d_a(4) E_field(:,2)*d_b(4) E_field(:,3)*d_c(4)];
    E_field_aligned = (inv(M)*deltaU')';
end

