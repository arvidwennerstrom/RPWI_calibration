clear all
% Position of Langmuir probes in spacecraft cooridnates [m].
% Rows for probe number [1-4] and columns for SC-axis [X Y Z].
% Source: "The Radio & Plasma Wave Investigation (RPWI)
% for the JUpiter ICy moons Explorer (JUICE), page 33)"
LP_positions = 1e-3*[2484 2895 4377; 1455 -3238 5303; 
                1455 -3238 -1849; -2768 2686 4432];

% Distance between probes in [X Y Z] directions. 4th column is norm of
% distance.
d_12 = LP_positions(1,:)-LP_positions(2,:); d_12 = [d_12 norm(d_12)];
d_13 = LP_positions(1,:)-LP_positions(3,:); d_13 = [d_13 norm(d_13)];
d_14 = LP_positions(1,:)-LP_positions(4,:); d_14 = [d_14 norm(d_14)];
d_23 = LP_positions(2,:)-LP_positions(3,:); d_23 = [d_23 norm(d_23)];
d_24 = LP_positions(2,:)-LP_positions(4,:); d_24 = [d_24 norm(d_24)];
d_34 = LP_positions(3,:)-LP_positions(4,:); d_34 = [d_34 norm(d_34)];
LP_delta_d = struct('d_12',d_12,'d_13',d_13,'d_14',d_14,'d_23',d_23,'d_24',d_24,'d_34',d_34);

save("C:\Users\1\Onedrive - KTH\MEX\IRF\Matlab\data\rpwi_data.mat","LP_positions","LP_delta_d")
% save("C:\Users\arvwe\Onedrive - KTH\MEX\IRF\Matlab\data\rpwi_data.mat","LP_positions","LP_delta_d")
% save("C:\Users\arvidwen\Onedrive - KTH\MEX\IRF\Matlab\rpwi_data.mat","LP_positions","delta_d")