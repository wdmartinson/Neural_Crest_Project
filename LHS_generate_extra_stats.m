function LHS_generate_extra_stats
% Program that allows you to collect perform eFAST sensitivity analysis on
% other summary statistics besides stream length, width, and nearest
% distances to cell

%% Make sure you are using the correct version of Python and g++
str = computer;
% desktop = true;
save_path = '/scratch/martinson/Documents/PhysiCell/NCC_puncta_model_output_data/';
if str(1) == 'M'
%     desktop = false;
    setenv('PATH', '/usr/local/bin:/usr/local/sbin:/usr/bin:/bin:/usr/sbin:/sbin');
    save_path = '/Volumes/easystore/NCC_puncta_model_output_data/';
end

%% Make vectors for the output data, so you can later analyse them:
% Get the total number of parameter regimes tested from the eFAST parameter
% matrix:
parameters.names = { 'fibronectin_puncta_wavelength',...
                     'filopodia_sensing_radius',...
                     'bias_towards_VM_direction',...
                     'average_time_until_next_filopodia_drop',...
                     'half_life_of_puncta_orientation_change',...
                     'y_length_of_cell_entrance_strip',...
                     'default_cell_radius',...
                     'default_cell_speed',...
                     'followers_cell_cell_repulsion_strength'};
parameters.values = readmatrix(strcat(save_path, 'LHS_Parameter_Values_Matrix.csv'));
number_of_parameter_regimes = size(parameters.values,1);

% Construct the data vectors:
output_average_cell_speed_in_stream = zeros(1, number_of_parameter_regimes);
output_average_FN_puncta_orientation = zeros(1, number_of_parameter_regimes);
output_average_number_of_cells_in_50_micron_ball = zeros(1, number_of_parameter_regimes);
output_average_PH_Stat_until_stream_all_connected = zeros(1, number_of_parameter_regimes);
output_average_gap_statistic = zeros(1, number_of_parameter_regimes);
output_average_gap_statistic_std_1_error = zeros(1, number_of_parameter_regimes);
output_average_leader_order_parameter = zeros(1, number_of_parameter_regimes);
output_average_follower_order_parameter = zeros(1, number_of_parameter_regimes);


if ~isfile(strcat(save_path, 'LHS_output_extra_stats.mat'))
%% Extract the relevant summary statistics from each parameter file and add them 
for j = 1:number_of_parameter_regimes
    fprintf(['Parameter set number ',num2str(j), ' of ',num2str(number_of_parameter_regimes),'...\n']);
    if isfile(strcat(save_path,num2str(j, '%04.f'),'.tar.gz'))
        % Unzip the .targz file:
        system(['tar -xf ',strcat(save_path,num2str(j, '%04.f'),'.tar.gz'), ' summary_statistics.mat']);
        load('summary_statistics.mat', 'Summary_Statistics');
        output_average_cell_speed_in_stream(j) = mean(Summary_Statistics.Average_Speed);
        output_average_FN_puncta_orientation(j) = mean(Summary_Statistics.Average_FN_Puncta_Orientations);
        output_average_number_of_cells_in_50_micron_ball(j) = mean(Summary_Statistics.Average_Number_of_Cells_within_50_um_Ball_per_Realization);
        output_average_PH_Stat_until_stream_all_connected(j) = mean(Summary_Statistics.Persistence_Homology_Largest_Radius_Before_All_Connected);
        output_average_gap_statistic(j) = mean(Summary_Statistics.Gap_Statistic);
        output_average_gap_statistic_std_1_error(j) = mean(Summary_Statistics.Gap_Statistic_Standard_1_Error);
        output_average_leader_order_parameter(j) = mean(Summary_Statistics.Order_Parameter_Leaders);
        output_average_follower_order_parameter(j) = mean(Summary_Statistics.Order_Parameter_Followers);     
        system('rm -rf summary_statistics.mat');
        continue;
    else
        error('You are missing at least one data set!'); 
    end % if you have already generated this data set
end % for j
output_average_PH_Stat_until_stream_all_connected(isinf(output_average_PH_Stat_until_stream_all_connected)) = 0;
output_average_follower_order_parameter(isnan(output_average_follower_order_parameter)) = 0;
output_average_number_of_cells_in_50_micron_ball(isnan(output_average_number_of_cells_in_50_micron_ball)) = 0;

%% Save the extra statistics:
save('LHS_output_extra_stats.mat', 'output_average_cell_speed_in_stream', 'output_average_FN_puncta_orientation', 'output_average_number_of_cells_in_50_micron_ball', 'output_average_PH_Stat_until_stream_all_connected', 'output_average_gap_statistic','output_average_gap_statistic_std_1_error','output_average_leader_order_parameter','output_average_follower_order_parameter');
system(['mv LHS_output_extra_stats.mat ', save_path]);
else
    load(strcat(save_path, 'LHS_output_extra_stats.mat'), 'output_average_cell_speed_in_stream', 'output_average_FN_puncta_orientation', 'output_average_number_of_cells_in_50_micron_ball', 'output_average_PH_Stat_until_stream_all_connected', 'output_average_gap_statistic','output_average_gap_statistic_std_1_error','output_average_leader_order_parameter','output_average_follower_order_parameter');
end % if the output file does not already exist
%% Determine PRCCs and their significance:
% List all of the .mat files currently in the directory (this should only 
% be the sobol indices):
load(strcat(save_path, 'PRCC_data.mat'), 'PRCCs', 'output_average_distance_to_nearest_cell');
output_average_distance_to_nearest_cell(isnan(output_average_distance_to_nearest_cell)) = 0;

alpha = 0.01; % Cut-off for statistical significance
[PRCCs.prcc_distance_to_nearest_cell, PRCCs.uncorrected_p_vals_distance_to_nearest_cell, PRCCs.significant_params_distance_to_nearest_cell]=PRCC(parameters.values,output_average_distance_to_nearest_cell,720,parameters.names,alpha);
[PRCCs.prcc_output_average_cell_speed_in_stream, PRCCs.uncorrected_p_vals_output_average_cell_speed_in_stream, PRCCs.significant_params_output_average_cell_speed_in_stream]=PRCC(parameters.values,output_average_cell_speed_in_stream,720,parameters.names,alpha);
[PRCCs.prcc_output_average_FN_puncta_orientation, PRCCs.uncorrected_p_vals_output_average_FN_puncta_orientation, PRCCs.significant_params_output_average_FN_puncta_orientation]=PRCC(parameters.values,output_average_FN_puncta_orientation,720,parameters.names,alpha);
[PRCCs.prcc_average_number_of_cells_in_50_micron_ball, PRCCs.uncorrected_p_vals_average_number_of_cells_in_50_micron_ball, PRCCs.significant_params_average_number_of_cells_in_50_micron_ball]=PRCC(parameters.values,output_average_number_of_cells_in_50_micron_ball,720,parameters.names,alpha);
[PRCCs.prcc_average_PH_Stat_until_stream_all_connected, PRCCs.uncorrected_p_vals_average_PH_Stat_until_stream_all_connected, PRCCs.significant_params_average_PH_Stat_until_stream_all_connected]=PRCC(parameters.values,output_average_PH_Stat_until_stream_all_connected,720,parameters.names,alpha);
[PRCCs.prcc_output_average_gap_statistic, PRCCs.uncorrected_p_vals_output_average_gap_statistic, PRCCs.significant_params_output_average_gap_statistic]=PRCC(parameters.values,output_average_gap_statistic,720,parameters.names,alpha);
[PRCCs.prcc_output_average_gap_statistic_std_1_error, PRCCs.uncorrected_p_vals_output_average_gap_statistic_std_1_error, PRCCs.significant_params_output_average_gap_statistic_std_1_error]=PRCC(parameters.values,output_average_gap_statistic_std_1_error,720,parameters.names,alpha);
[PRCCs.prcc_output_average_leader_order_parameter, PRCCs.uncorrected_p_vals_output_average_leader_order_parameter, PRCCs.significant_params_output_average_leader_order_parameter]=PRCC(parameters.values,output_average_leader_order_parameter,720,parameters.names,alpha);
[PRCCs.prcc_output_average_follower_order_parameter, PRCCs.uncorrected_p_vals_output_average_follower_order_parameter, PRCCs.significant_params_output_average_follower_order_parameter]=PRCC(parameters.values,output_average_follower_order_parameter,720,parameters.names,alpha);

save('PRCCs_with_extra_stats.mat', 'PRCCs');
system(['mv ./PRCCs_with_extra_stats.mat ', save_path]);
end