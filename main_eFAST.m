% LHS Global Sensitivity Analysis Code for NCC Puncta model
%% Make sure you are using the correct version of Python and g++
str = computer;
desktop = true;
save_path = './NCC_puncta_model_output_data_eFAST/';
if str(1) == 'M'
    desktop = false;
    setenv('PATH', '/usr/local/bin:/usr/local/sbin:/usr/bin:/bin:/usr/sbin:/sbin');
    save_path = './NCC_puncta_model_output_data_eFAST/';
end
%% Create cell with all of the parameter values you wish to change
parameters.names =  {'filopodia_sensing_radius',...
                     'bias_towards_VM_direction',...
                     'followers_cell_cell_repulsion_strength',...
                     'dummy_parameter'};
                 
%% Create Latin Hypercube Sampling Matrix, adapting open-source code from Marino et al. (2008, J. Theor. Biol.)
number_of_parameters = length(parameters.names);

% Write the parameter value matrix as a MAT file for safe keeping:
if ~isfile(strcat(save_path, 'eFAST_sample_matrix.mat'))
    system(['python3 eFAST_generate_samples.py ',save_path, ' && cp eFAST_sample_matrix.mat ', save_path]);
    load(strcat(save_path, 'eFAST_sample_matrix.mat'));
else
    load(strcat(save_path, 'eFAST_sample_matrix.mat'));
end
% eFAST_sample_matrix is a 3D array. Reformat to make it 2D so that it's
% easier to extract output:
parameters.values = zeros(size(eFAST_sample_matrix,1)*size(eFAST_sample_matrix,3), size(eFAST_sample_matrix,2));
for z = 1:size(eFAST_sample_matrix,3)
    parameters.values((z-1)*size(eFAST_sample_matrix,1)+1:z*size(eFAST_sample_matrix,1), :) = eFAST_sample_matrix(:,:,z);
end
number_of_parameter_regimes = size(parameters.values, 1);
number_of_realizations_per_parameter_regime = 20; 
maximum_realization_time_per_realization = 1.5; % min
est_runtime_hours = number_of_parameter_regimes*(number_of_realizations_per_parameter_regime*maximum_realization_time_per_realization+1)/60; % Add on an extra minute for producing the video files and generating the summary statistics
est_runtime_days = est_runtime_hours/24;

fprintf(strcat('The estimated maximum amount of time that this parameter sweep will take to complete is \t', num2str(est_runtime_hours),' hours or \t', num2str(est_runtime_days), ' days.\n'));
input('Do you still wish to continue (press Enter if yes, Ctrl+C if no)?');

Values_for_this_parameter_set = cell(1,number_of_parameters-1);

output_average_stream_length = zeros(number_of_parameter_regimes, 1); % Output vector for estimated stream length
output_average_stream_width = zeros(number_of_parameter_regimes, 1); % Output vector for estimated stream width
output_average_max_cell_x_location = zeros(number_of_parameter_regimes, 1); % Output vector for estimated stream length with no FN orientation
output_average_max_cell_y_width = zeros(number_of_parameter_regimes, 1); % Output vector for estimated stream width with no FN orientation
output_average_distance_to_nearest_cell = zeros(number_of_parameter_regimes, 1); % Output vector for nearest distance to a cell
%% Loop over each parameter value, edit the XML config file, and run for multiple realizations:
for j = 1:number_of_parameter_regimes
    for col = 1:number_of_parameters-1 % don't count the dummy parameter
        Values_for_this_parameter_set{col} = parameters.values(j, col);
    end % for col
    if isfile(strcat(save_path,num2str(j, '%04.f'),'.tar.gz'))
        % Unzip the .targz file:
        system(['tar -xf ',strcat(save_path,num2str(j, '%04.f'),'.tar.gz'), ' summary_statistics.mat']);
        load('summary_statistics.mat');
        output_average_stream_length(j) = mean(Summary_Statistics.Est_Stream_Length);
        output_average_stream_width(j) = mean(Summary_Statistics.Est_Stream_Width);
        output_average_max_cell_x_location(j) = mean(Summary_Statistics.Max_X_Length_of_Stream);
        output_average_max_cell_y_width(j) = mean(Summary_Statistics.Max_Y_Length_of_Stream);
        output_average_distance_to_nearest_cell(j) = mean(Summary_Statistics.Average_Min_Distance_per_Realization);
        system('rm -rf summary_statistics.mat');
        continue;
    end % if you have already generated this data set

    Parameters = {parameters.names{1:number_of_parameters-1}; Values_for_this_parameter_set{1:number_of_parameters-1}}; % don't count the dummy parameter
    py_file = ' ./edit_PhysiCell_settings_XML.py';
    cmd = strcat('python3 ', py_file, sprintf(' %s %d', Parameters{:}));  
    % Use a python function to edit the config file
    system(cmd);       
    fprintf(['Parameter set number ',num2str(j), ' of ',num2str(number_of_parameter_regimes),'...\n']);
    %% Compute all of the summary statistics you require:                         
    Summary_Statistics = NCC_Model_Puncta_Multiple_Realizations_eFAST(j, number_of_realizations_per_parameter_regime);
    output_average_stream_length(j) = mean(Summary_Statistics.Est_Stream_Length);
    output_average_stream_width(j) = mean(Summary_Statistics.Est_Stream_Width);
    output_average_max_cell_x_location(j) = mean(Summary_Statistics.Max_X_Length_of_Stream);
    output_average_max_cell_y_width(j) = mean(Summary_Statistics.Max_Y_Length_of_Stream);
    output_average_distance_to_nearest_cell(j) = mean(Summary_Statistics.Average_Min_Distance_per_Realization);
end % for j
output_average_distance_to_nearest_cell(isnan(output_average_distance_to_nearest_cell)) = 0;
% Reformat the output to 3D arrays to match the size of
% eFAST_sample_matrix:
output_average_stream_length = reshape(output_average_stream_length, size(eFAST_sample_matrix,1), 1, size(eFAST_sample_matrix,3) );
output_average_stream_width = reshape(output_average_stream_width, size(eFAST_sample_matrix,1), 1, size(eFAST_sample_matrix,3) );
output_average_max_cell_x_location = reshape(output_average_max_cell_x_location, size(eFAST_sample_matrix,1), 1, size(eFAST_sample_matrix,3) );
output_average_max_cell_y_width = reshape(output_average_max_cell_y_width, size(eFAST_sample_matrix,1), 1, size(eFAST_sample_matrix,3) );
output_average_distance_to_nearest_cell = reshape(output_average_distance_to_nearest_cell, size(eFAST_sample_matrix,1), 1, size(eFAST_sample_matrix,3) );
save('eFAST_output_stats.mat', 'output_average_stream_length', 'output_average_stream_width', 'output_average_max_cell_x_location', 'output_average_max_cell_y_width', 'output_average_distance_to_nearest_cell');
%% Compute Sobol Indicies
system('python3 eFAST_analyze_results.py');
system(['mv eFAST_output_stats.mat ', save_path, ' && rm -rf eFAST_sample_matrix.mat']);
%% Determine significance using 2-sided t-test:
% List all of the .mat files currently in the directory (this should only 
% be the sobol indices):
listing = dir('./*.mat');
N = length(listing);
alpha = 0.01; % cutoff value for statistical significance
for i = 1:N
    load(listing(i).name);
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