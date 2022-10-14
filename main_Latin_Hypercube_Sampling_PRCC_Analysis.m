% LHS Global Sensitivity Analysis Code for NCC Puncta model
%% Make sure you are using the correct version of Python and g++
str = computer;
desktop = true;
if str(1) == 'M'
    desktop = false;
    setenv('PATH', '/usr/local/bin:/usr/local/sbin:/usr/bin:/bin:/usr/sbin:/sbin');
end
%% Create cell with all of the parameter values you wish to change
parameters.names = { 'fibronectin_puncta_wavelength',...
                     'filopodia_sensing_radius',...
                     'bias_towards_VM_direction',...
                     'average_time_until_next_filopodia_drop',...
                     'half_life_of_puncta_orientation_change',...
                     'y_length_of_cell_entrance_strip',...
                     'default_cell_radius',...
                     'default_cell_speed',...
                     'followers_cell_cell_repulsion_strength'};
                 
%% Create Latin Hypercube Sampling Matrix, adapting open-source code from Marino et al. (2008, J. Theor. Biol.)
number_of_parameters = length(parameters.names);
number_of_parameter_regimes = 1000; % Sample 1000 different parameter values
parameters.values = zeros(number_of_parameter_regimes, number_of_parameters);
% PARAMETER BASELINE VALUES
lambda_FN = 30;
R_filo = 40;
rho = 0.5;
T_ave = 30;
T_half = 60;
l_entr = 50;
R_cell = 10;
s_FN = 0.5;
c_i = 1;

% Use Latin Hypercube Sampling to Select the Parameter Values:
rng(1); % for reproducabiltiy

% parameters.values(:,1)=LHS_Call(20, lambda_FN, 75, number_of_parameter_regimes,'unif'); 
% parameters.values(:,2)=LHS_Call(20, R_filo, 75, number_of_parameter_regimes,'unif'); 
% parameters.values(:,3)=LHS_Call(0.25, rho, 1, number_of_parameter_regimes,'unif'); 
% parameters.values(:,4)=LHS_Call(15, T_ave, 90, number_of_parameter_regimes,'unif'); 
% parameters.values(:,5)=LHS_Call(5, T_half, 90, number_of_parameter_regimes,'unif');
% parameters.values(:,6)=LHS_Call(30, l_entr, 500, number_of_parameter_regimes,'unif');
% parameters.values(:,7)=LHS_Call(7.5, R_cell, 17.5, number_of_parameter_regimes,'unif');
% parameters.values(:,8)=LHS_Call(0.5, s_FN, 1.25, number_of_parameter_regimes,'unif');
% parameters.values(:,9)=LHS_Call(1e-5, c_i, 5, number_of_parameter_regimes,'unif');
% parameters.values(:,10)= lhsdesign(number_of_parameter_regimes,1);

parameters.values = lhsdesign(number_of_parameter_regimes, number_of_parameters); % Generate LHS matrix from Unif([0,1]) distribution
parameters.values(:,1:2) = 20 + 55*parameters.values(:,1:2); % % lambda_FN = [20, 75], R_filo = [20, 75]
parameters.values(:,3) = 0.25 + 0.75*parameters.values(:,3); % rho = [0.25, 1]
parameters.values(:,4) = 15 + 75*parameters.values(:,4); % T_ave = [15, 90]
parameters.values(:,5) = 2 + 88*parameters.values(:,5); % T_half = [2, 90]
parameters.values(:,6) = 30 + 470*parameters.values(:,6); % l_entr = [30, 500]
parameters.values(:,7) = 7.5 + 10*parameters.values(:,7); % R_cell = [7.5, 17.5]
parameters.values(:,8) = 0.5 + 0.75*parameters.values(:,8); % s_FN = [0.5, 1.25]
parameters.values(:,9) = 1e-5 + (5-1e-5)*parameters.values(:,9); % c_i = [1e-5, 5]

% Write the LHS matrix as a CSV file for safe keeping:
if ~isfile('./NCC_puncta_model_LHS_data/LHS_Parameter_Values_Matrix.csv')
    writematrix(parameters.values, './NCC_puncta_model_LHS_data/LHS_Parameter_Values_Matrix.csv');
else
    parameters.values = readmatrix('./NCC_puncta_model_LHS_data/LHS_Parameter_Values_Matrix.csv');
end
number_of_realizations_per_parameter_regime = 20;
est_maximum_realization_time_per_realization = 1.5; % min
est_runtime_hours = number_of_parameter_regimes*(number_of_realizations_per_parameter_regime*est_maximum_realization_time_per_realization+1)/60; % Add on an extra minute for producing the video files and generating the summary statistics
est_runtime_days = est_runtime_hours/24;

fprintf(strcat('The estimated maximum amount of time that this parameter sweep will take to complete is \t', num2str(est_runtime_hours),' hours or \t', num2str(est_runtime_days), ' days.\n'));
input('Do you still wish to continue (press Enter if yes, Ctrl+C if no)?');

Values_for_this_parameter_set = cell(1,number_of_parameters);

output_average_stream_length = zeros(1, number_of_parameter_regimes); % Output vector for estimated stream length
output_average_stream_width = zeros(1, number_of_parameter_regimes); % Output vector for estimated stream width
output_average_max_cell_x_location = zeros(1, number_of_parameter_regimes); % Output vector for estimated stream length with no FN orientation
output_average_max_cell_y_width = zeros(1, number_of_parameter_regimes); % Output vector for estimated stream width with no FN orientation
output_average_distance_to_nearest_cell = zeros(1, number_of_parameter_regimes); % Output vector for nearest distance to a cell
%% Loop over each parameter value, edit the XML config file, and run for multiple realizations:
for j = 1:number_of_parameter_regimes
    for col = 1:number_of_parameters
        Values_for_this_parameter_set{col} = parameters.values(j, col);
    end % for col
    if isfile(strcat('./NCC_puncta_model_LHS_data/',num2str(j, '%04.f'),'.tar.gz')) 
        % Unzip the .targz file:
        system(['tar -xf ',strcat('./NCC_puncta_model_LHS_data/',num2str(j, '%04.f'),'.tar.gz'), ' summary_statistics.mat']);
        load('summary_statistics.mat');
        output_average_stream_length(j) = mean(Summary_Statistics.Est_Stream_Length);
        output_average_stream_width(j) = mean(Summary_Statistics.Est_Stream_Width);
        output_average_max_cell_x_location(j) = mean(Summary_Statistics.Max_X_Length_of_Stream);
        output_average_max_cell_y_width(j) = mean(Summary_Statistics.Max_Y_Length_of_Stream);
        output_average_distance_to_nearest_cell(j) = mean(Summary_Statistics.Average_Min_Distance_per_Realization);
        system('rm -rf summary_statistics.mat');
        continue;
    end % if you have already generated this data set

    Parameters = {parameters.names{1:number_of_parameters}; Values_for_this_parameter_set{1:number_of_parameters}};
    py_file = ' ./edit_PhysiCell_settings_XML.py';
    cmd = strcat('python3 ', py_file, sprintf(' %s %d', Parameters{:}));  
    % Use a python function to edit the config file
    system(cmd);       
    fprintf(['Parameter set number ',num2str(j), ' of ',num2str(number_of_parameter_regimes),'...\n']);
    %% Compute all of the summary statistics you require:                         
    Summary_Statistics = NCC_Model_Puncta_Multiple_Realizations(j, number_of_realizations_per_parameter_regime);
    output_average_stream_length(j) = mean(Summary_Statistics.Est_Stream_Length);
    output_average_stream_width(j) = mean(Summary_Statistics.Est_Stream_Width);
    output_average_max_cell_x_location(j) = mean(Summary_Statistics.Max_X_Length_of_Stream);
    output_average_max_cell_y_width(j) = mean(Summary_Statistics.Max_Y_Length_of_Stream);
    output_average_distance_to_nearest_cell(j) = mean(Summary_Statistics.Average_Min_Distance_per_Realization);
end % for j

%% Compute PRCCs using the PRCC.m script:
[PRCCs.prcc_est_stream_length, PRCCs.uncorrected_p_vals_est_stream_length, PRCCs.significant_params_est_stream_length]=PRCC(parameters.values,output_average_stream_length,750,parameters.names,0.01);
[PRCCs.prcc_est_stream_width, PRCCs.uncorrected_p_vals_est_stream_width, PRCCs.significant_params_est_stream_width]=PRCC(parameters.values,output_average_stream_width,750,parameters.names,0.01);
[PRCCs.prcc_max_cell_x_location, PRCCs.uncorrected_p_vals_max_cell_x_location, PRCCs.significant_params_max_cell_x_location]=PRCC(parameters.values,output_average_max_cell_x_location,750,parameters.names,0.01);
[PRCCs.prcc_max_cell_y_width, PRCCs.uncorrected_p_vals_max_cell_y_width, PRCCs.significant_params_max_cell_y_width]=PRCC(parameters.values,output_average_max_cell_y_width,750,parameters.names,0.01);
[PRCCs.prcc_distance_to_nearest_cell, PRCCs.uncorrected_p_vals_distance_to_nearest_cell, PRCCs.significant_params_distance_to_nearest_cell]=PRCC(parameters.values,output_average_distance_to_nearest_cell,750,parameters.names,0.01);
save('./NCC_puncta_model_LHS_data/PRCC_data.mat', 'PRCCs', 'output_average_stream_length', 'output_average_stream_width', 'output_average_max_cell_x_location', 'output_average_max_cell_y_width', 'output_average_distance_to_nearest_cell');
