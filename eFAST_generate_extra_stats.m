function eFAST_generate_extra_stats
% Program that allows you to collect perform eFAST sensitivity analysis on
% other summary statistics besides stream length, width, and nearest
% distances to cell

%% Make sure you are using the correct version of Python and g++
str = computer;
% desktop = true;
save_path = '/scratch/martinson/Documents/PhysiCell/NCC_puncta_model_output_data_eFAST/';
if str(1) == 'M'
%     desktop = false;
    setenv('PATH', '/usr/local/bin:/usr/local/sbin:/usr/bin:/bin:/usr/sbin:/sbin');
    save_path = '/Volumes/easystore/NCC_puncta_model_output_data_eFAST/';
end

%% Make vectors for the output data, so you can later analyse them:
% Get the total number of parameter regimes tested from the eFAST parameter
% matrix:
load(strcat(save_path, 'eFAST_sample_matrix.mat'), 'eFAST_sample_matrix');
system(['cp ',strcat(save_path, 'eFAST_sample_matrix.mat'), ' ./']);
number_of_parameter_regimes = size(eFAST_sample_matrix,1)*size(eFAST_sample_matrix,3);
number_of_parameters = size(eFAST_sample_matrix,2);

% Construct the data vectors:
output_average_cell_speed_in_stream = zeros(number_of_parameter_regimes, 1);
output_average_FN_puncta_orientation = zeros(number_of_parameter_regimes, 1);
output_average_number_of_cells_in_50_micron_ball = zeros(number_of_parameter_regimes, 1);
output_average_PH_Stat_until_stream_all_connected = zeros(number_of_parameter_regimes, 1);
output_average_gap_statistic = zeros(number_of_parameter_regimes, 1);
output_average_gap_statistic_std_1_error = zeros(number_of_parameter_regimes, 1);
output_average_leader_order_parameter = zeros(number_of_parameter_regimes, 1);
output_average_follower_order_parameter = zeros(number_of_parameter_regimes, 1);


if ~isfile(strcat(save_path, 'eFAST_output_extra_stats.mat'))
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

% Reformat the output to 3D arrays to match the size of
% eFAST_sample_matrix:
output_average_cell_speed_in_stream = reshape(output_average_cell_speed_in_stream, size(eFAST_sample_matrix,1), 1, size(eFAST_sample_matrix,3) );
output_average_FN_puncta_orientation = reshape(output_average_FN_puncta_orientation, size(eFAST_sample_matrix,1), 1, size(eFAST_sample_matrix,3) );
output_average_number_of_cells_in_50_micron_ball = reshape(output_average_number_of_cells_in_50_micron_ball, size(eFAST_sample_matrix,1), 1, size(eFAST_sample_matrix,3) );
output_average_PH_Stat_until_stream_all_connected = reshape(output_average_PH_Stat_until_stream_all_connected, size(eFAST_sample_matrix,1), 1, size(eFAST_sample_matrix,3) );
output_average_gap_statistic = reshape(output_average_gap_statistic, size(eFAST_sample_matrix,1), 1, size(eFAST_sample_matrix,3) );
output_average_gap_statistic_std_1_error = reshape(output_average_gap_statistic_std_1_error, size(eFAST_sample_matrix,1), 1, size(eFAST_sample_matrix,3) );
output_average_leader_order_parameter = reshape(output_average_leader_order_parameter, size(eFAST_sample_matrix,1), 1, size(eFAST_sample_matrix,3) );
output_average_follower_order_parameter = reshape(output_average_follower_order_parameter, size(eFAST_sample_matrix,1), 1, size(eFAST_sample_matrix,3) );

%% Save the extra statistics and compute Sobol indices:
save('eFAST_output_extra_stats.mat', 'output_average_cell_speed_in_stream', 'output_average_FN_puncta_orientation', 'output_average_number_of_cells_in_50_micron_ball', 'output_average_PH_Stat_until_stream_all_connected', 'output_average_gap_statistic','output_average_gap_statistic_std_1_error','output_average_leader_order_parameter','output_average_follower_order_parameter');
system('python3 eFAST_analyze_results.py');
system(['mv eFAST_output_extra_stats.mat ', save_path, ' && rm -rf eFAST_output_extra_stats.mat && rm -rf eFAST_sample_matrix.mat']);
else
system(['cp ',strcat(save_path, 'eFAST_output_extra_stats.mat'), ' ./ && python3 eFAST_analyze_results.py && rm -rf eFAST_output_extra_stats.mat && rm -rf eFAST_sample_matrix.mat']);
end % if the output file does not already exist
%% Determine significance using 2-sided t-test:
% List all of the .mat files currently in the directory (this should only 
% be the sobol indices):
listing = dir('./*.mat');
N = length(listing);
alpha = 0.01; % cutoff value for statistical significance
for i = 1:N
    load(listing(i).name, 'sobol_indices');
    M = length(sobol_indices);
    Si_vals = zeros(M, number_of_parameters);
    St_vals = zeros(M, number_of_parameters);
    for jj = 1:M
        Si_vals(jj, :) = sobol_indices{1,jj}.S1;
        St_vals(jj, :) = sobol_indices{1,jj}.ST;
    end % for jj
    Si_are_you_signficant = zeros(number_of_parameters-1, 1);
    Si_pvals_significance = zeros(number_of_parameters-1, 1);
    St_are_you_signficant = zeros(number_of_parameters-1, 1);
    St_pvals_significance = zeros(number_of_parameters-1, 1);
    for pp = 1:number_of_parameters-1
        % CV Test ? check that num_resamples is sufficient 
        %   (okay for Si, St >> 0)
        CV_Si = std(Si_vals, 0, 1)./mean(Si_vals, 1);
        CV_St = std(St_vals, 0, 1)./mean(St_vals, 1);
        % 2-sided t-test for statistical significance:
        [Si_are_you_signficant(pp), Si_pvals_significance(pp)] = ttest2(Si_vals(:,pp), Si_vals(:, end), alpha , 'right', 'unequal');
        [St_are_you_signficant(pp), St_pvals_significance(pp)] = ttest2(St_vals(:,pp), St_vals(:, end), alpha , 'right', 'unequal');
    end
    save(listing(i).name, 'sobol_indices', 'Si_vals', 'St_vals', 'Si_are_you_signficant', 'Si_pvals_significance', 'St_are_you_signficant', 'St_pvals_significance', 'CV_Si', 'CV_St');
end % for i
system(['mv ./*.mat ', save_path]);
end