function eFAST_compute_PRCCs
% Compute PRCCs of eFAST results, so that you can compare with your LHS
% analysis.

%% Make sure you are using the correct version of Python and g++
str = computer;
% desktop = true;
save_path = '/scratch/martinson/Documents/PhysiCell/NCC_puncta_model_output_data_eFAST/';
if str(1) == 'M'
%     desktop = false;
    setenv('PATH', '/usr/local/bin:/usr/local/sbin:/usr/bin:/bin:/usr/sbin:/sbin');
    save_path = '/Volumes/easystore/NCC_puncta_model_output_data_eFAST/';
end

parameters.names =  {'filopodia_sensing_radius',...
                     'bias_towards_VM_direction',...
                     'followers_cell_cell_repulsion_strength',...
                     'dummy_parameter'};

load([save_path,'eFAST_sample_matrix.mat'], 'eFAST_sample_matrix');
load([save_path,'eFAST_output_stats.mat'], 'output_average_stream_length', 'output_average_stream_width', 'output_average_max_cell_x_location', 'output_average_max_cell_y_width' , 'output_average_distance_to_nearest_cell');
load([save_path,'eFAST_output_extra_stats.mat'], 'output_average_cell_speed_in_stream', 'output_average_FN_puncta_orientation', 'output_average_number_of_cells_in_50_micron_ball', 'output_average_PH_Stat_until_stream_all_connected', 'output_average_gap_statistic','output_average_gap_statistic_std_1_error','output_average_leader_order_parameter','output_average_follower_order_parameter');

% Reshape all matrices to vectors
parameters.values = zeros(size(eFAST_sample_matrix,1)*size(eFAST_sample_matrix,3), size(eFAST_sample_matrix,2));
for z = 1:size(eFAST_sample_matrix,3)
    parameters.values((z-1)*size(eFAST_sample_matrix,1)+1:z*size(eFAST_sample_matrix,1), :) = eFAST_sample_matrix(:,:,z);
end
output_average_stream_length = reshape(output_average_stream_length, size(eFAST_sample_matrix,1)*size(eFAST_sample_matrix,3), 1 )';
output_average_stream_width = reshape(output_average_stream_width, size(eFAST_sample_matrix,1)*size(eFAST_sample_matrix,3), 1 )';
output_average_max_cell_x_location = reshape(output_average_max_cell_x_location, size(eFAST_sample_matrix,1)*size(eFAST_sample_matrix,3), 1 )';
output_average_max_cell_y_width = reshape(output_average_max_cell_y_width, size(eFAST_sample_matrix,1)*size(eFAST_sample_matrix,3), 1 )';
output_average_distance_to_nearest_cell = reshape(output_average_distance_to_nearest_cell, size(eFAST_sample_matrix,1)*size(eFAST_sample_matrix,3), 1 )';
output_average_cell_speed_in_stream = reshape(output_average_cell_speed_in_stream, size(eFAST_sample_matrix,1)*size(eFAST_sample_matrix,3), 1 )';
output_average_FN_puncta_orientation = reshape(output_average_FN_puncta_orientation, size(eFAST_sample_matrix,1)*size(eFAST_sample_matrix,3), 1 )';
output_average_number_of_cells_in_50_micron_ball = reshape(output_average_number_of_cells_in_50_micron_ball, size(eFAST_sample_matrix,1)*size(eFAST_sample_matrix,3), 1 )';
output_average_PH_Stat_until_stream_all_connected = reshape(output_average_PH_Stat_until_stream_all_connected, size(eFAST_sample_matrix,1)*size(eFAST_sample_matrix,3), 1 )';
output_average_gap_statistic = reshape(output_average_gap_statistic, size(eFAST_sample_matrix,1)*size(eFAST_sample_matrix,3), 1 )';
output_average_gap_statistic_std_1_error = reshape(output_average_gap_statistic_std_1_error, size(eFAST_sample_matrix,1)*size(eFAST_sample_matrix,3), 1 )';
output_average_leader_order_parameter = reshape(output_average_leader_order_parameter, size(eFAST_sample_matrix,1)*size(eFAST_sample_matrix,3), 1 )';
output_average_follower_order_parameter = reshape(output_average_follower_order_parameter, size(eFAST_sample_matrix,1)*size(eFAST_sample_matrix,3), 1 )';


alpha = 0.01; % Cut-off for statistical significance
[PRCCs.prcc_est_stream_length, PRCCs.uncorrected_p_vals_est_stream_length, PRCCs.significant_params_est_stream_length]=PRCC(parameters.values,output_average_stream_length,720,parameters.names,alpha);
[PRCCs.prcc_est_stream_width, PRCCs.uncorrected_p_vals_est_stream_width, PRCCs.significant_params_est_stream_width]=PRCC(parameters.values,output_average_stream_width,720,parameters.names,alpha);
[PRCCs.prcc_max_cell_x_location, PRCCs.uncorrected_p_vals_max_cell_x_location, PRCCs.significant_params_max_cell_x_location]=PRCC(parameters.values,output_average_max_cell_x_location,720,parameters.names,alpha);
[PRCCs.prcc_max_cell_y_width, PRCCs.uncorrected_p_vals_max_cell_y_width, PRCCs.significant_params_max_cell_y_width]=PRCC(parameters.values,output_average_max_cell_y_width,720,parameters.names,alpha);
[PRCCs.prcc_distance_to_nearest_cell, PRCCs.uncorrected_p_vals_distance_to_nearest_cell, PRCCs.significant_params_distance_to_nearest_cell]=PRCC(parameters.values,output_average_distance_to_nearest_cell,720,parameters.names,alpha);
[PRCCs.prcc_distance_to_nearest_cell, PRCCs.uncorrected_p_vals_distance_to_nearest_cell, PRCCs.significant_params_distance_to_nearest_cell]=PRCC(parameters.values,output_average_distance_to_nearest_cell,720,parameters.names,alpha);
[PRCCs.prcc_output_average_cell_speed_in_stream, PRCCs.uncorrected_p_vals_output_average_cell_speed_in_stream, PRCCs.significant_params_output_average_cell_speed_in_stream]=PRCC(parameters.values,output_average_cell_speed_in_stream,720,parameters.names,alpha);
[PRCCs.prcc_output_average_FN_puncta_orientation, PRCCs.uncorrected_p_vals_output_average_FN_puncta_orientation, PRCCs.significant_params_output_average_FN_puncta_orientation]=PRCC(parameters.values,output_average_FN_puncta_orientation,720,parameters.names,alpha);
[PRCCs.prcc_average_number_of_cells_in_50_micron_ball, PRCCs.uncorrected_p_vals_average_number_of_cells_in_50_micron_ball, PRCCs.significant_params_average_number_of_cells_in_50_micron_ball]=PRCC(parameters.values,output_average_number_of_cells_in_50_micron_ball,720,parameters.names,alpha);
[PRCCs.prcc_average_PH_Stat_until_stream_all_connected, PRCCs.uncorrected_p_vals_average_PH_Stat_until_stream_all_connected, PRCCs.significant_params_average_PH_Stat_until_stream_all_connected]=PRCC(parameters.values,output_average_PH_Stat_until_stream_all_connected,720,parameters.names,alpha);
[PRCCs.prcc_output_average_gap_statistic, PRCCs.uncorrected_p_vals_output_average_gap_statistic, PRCCs.significant_params_output_average_gap_statistic]=PRCC(parameters.values,output_average_gap_statistic,720,parameters.names,alpha);
[PRCCs.prcc_output_average_gap_statistic_std_1_error, PRCCs.uncorrected_p_vals_output_average_gap_statistic_std_1_error, PRCCs.significant_params_output_average_gap_statistic_std_1_error]=PRCC(parameters.values,output_average_gap_statistic_std_1_error,720,parameters.names,alpha);
[PRCCs.prcc_output_average_leader_order_parameter, PRCCs.uncorrected_p_vals_output_average_leader_order_parameter, PRCCs.significant_params_output_average_leader_order_parameter]=PRCC(parameters.values,output_average_leader_order_parameter,720,parameters.names,alpha);
[PRCCs.prcc_output_average_follower_order_parameter, PRCCs.uncorrected_p_vals_output_average_follower_order_parameter, PRCCs.significant_params_output_average_follower_order_parameter]=PRCC(parameters.values,output_average_follower_order_parameter,720,parameters.names,alpha);

save([save_path, 'PRCCs_from_eFAST.mat'], 'PRCCs');
end