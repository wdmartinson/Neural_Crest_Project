% Consistency analysis for NCC ABM
%% Make sure you are using the correct version of Python and g++
str = computer;
desktop = true;
save_path = './NCC_consistency_analysis/';
if str(1) == 'M'
    desktop = false;
    setenv('PATH', '/usr/local/bin:/usr/local/sbin:/usr/bin:/bin:/usr/sbin:/sbin');
    save_path = './NCC_consistency_analysis/';
end
close all;
%% Set parameter values:
parameters.names = { 'fibronectin_puncta_wavelength',...
                     'filopodia_sensing_radius',...
                     'bias_towards_VM_direction',...
                     'average_time_until_next_filopodia_drop',...
                     'half_life_of_puncta_orientation_change',...
                     'y_length_of_cell_entrance_strip',...
                     'default_cell_radius',...
                     'default_cell_speed',...
                     'followers_cell_cell_repulsion_strength'};
% PARAMETER VALUES for NCC model (known to give successful stream
%   formation)
lambda_FN = 35;
R_filo = 50;
rho = 0.5;
T_ave = 30;
T_half = 5;
l_entr = 50;
R_cell = 10;
s_FN = 0.5;
c_i = 1;

parameters.values = cell(1, length(parameters.names));
parameters.values{1} = lambda_FN;
parameters.values{2} = R_filo;
parameters.values{3} = rho;
parameters.values{4} = T_ave;
parameters.values{5} = T_half;
parameters.values{6} = l_entr;
parameters.values{7} = R_cell;
parameters.values{8} = s_FN;
parameters.values{9} = c_i;

%% Set up hyperparameters
number_of_parameters = length(parameters.names);
number_of_distributions = 20;
number_of_realizations_per_distribution = [1, 20, 50, 100, 200];
number_of_groups = length(number_of_realizations_per_distribution);
maximum_realization_time_per_realization = 1.5; % min
est_runtime_hours = number_of_distributions*sum(number_of_realizations_per_distribution)*maximum_realization_time_per_realization/60;
est_runtime_days = est_runtime_hours/24;

fprintf(strcat('The estimated maximum amount of time that this consistency analysis will take to complete is \t', num2str(est_runtime_hours),' hours or \t', num2str(est_runtime_days), ' days.\n'));
input('Do you still wish to continue (press Enter if yes, Ctrl+C if no)?');

% Create 3D array for storing the relevant values after each realization:
output_average_stream_length = zeros(max(number_of_realizations_per_distribution), number_of_distributions, number_of_groups); % Output array for estimated stream length
output_average_stream_width = zeros(max(number_of_realizations_per_distribution), number_of_distributions, number_of_groups); % Output array for estimated stream width
output_average_max_cell_x_location = zeros(max(number_of_realizations_per_distribution), number_of_distributions, number_of_groups); % Output array for estimated stream length with no FN orientation
output_average_max_cell_y_width = zeros(max(number_of_realizations_per_distribution), number_of_distributions, number_of_groups); % Output array for estimated stream width with no FN orientation
output_average_distance_to_nearest_cell = zeros(max(number_of_realizations_per_distribution), number_of_distributions, number_of_groups); % Output array for nearest distance to a cell

% Also check on the cumulative means:
cum_mean_average_stream_length = zeros(1,number_of_distributions*sum(number_of_realizations_per_distribution));
cum_mean_average_stream_width = zeros(1,number_of_distributions*sum(number_of_realizations_per_distribution));
cum_mean_average_max_cell_x_location = zeros(1,number_of_distributions*sum(number_of_realizations_per_distribution));
cum_mean_average_max_cell_y_width = zeros(1,number_of_distributions*sum(number_of_realizations_per_distribution));
cum_mean_average_distance_to_nearest_cell = zeros(1,number_of_distributions*sum(number_of_realizations_per_distribution));

% Use a python function to edit the config file
Parameters = {parameters.names{1:number_of_parameters}; parameters.values{1:number_of_parameters}};
py_file = ' ./edit_PhysiCell_settings_XML.py';
cmd = strcat('python3 ', py_file, sprintf(' %s %d', Parameters{:}));  
system(cmd);

%% Collect Data For Consistency Analysis
for G = 1:number_of_groups
    for D = 1:number_of_distributions
        video_title = strcat(num2str(number_of_realizations_per_distribution(G)),'_',num2str(D));
        if isfile(strcat(save_path,video_title,'.tar.gz'))  
            %  Unzip the .tar.gz file:
            system(['tar -xf ',strcat(save_path,video_title,'.tar.gz'), ' ',video_title,'_summary_statistics.mat']);
            load([video_title,'_summary_statistics.mat']);
            output_average_stream_length(1:number_of_realizations_per_distribution(G), D, G) = Summary_Statistics.Est_Stream_Length;
            output_average_stream_width(1:number_of_realizations_per_distribution(G), D, G) = Summary_Statistics.Est_Stream_Width;
            output_average_max_cell_x_location(1:number_of_realizations_per_distribution(G), D, G) = Summary_Statistics.Max_X_Length_of_Stream;
            output_average_max_cell_y_width(1:number_of_realizations_per_distribution(G), D, G) = Summary_Statistics.Max_Y_Length_of_Stream;
            output_average_distance_to_nearest_cell(1:number_of_realizations_per_distribution(G), D, G) = Summary_Statistics.Average_Min_Distance_per_Realization;
            
            for rr = 1:number_of_realizations_per_distribution(G)
                idx = number_of_distributions*sum(number_of_realizations_per_distribution(1:G-1)) + (D-1)*number_of_realizations_per_distribution(G) + rr;
                cum_mean_average_stream_length(idx) = ( sum(sum(sum(output_average_stream_length(1:number_of_realizations_per_distribution(G), :, 1:G-1)))) + sum(sum(sum(output_average_stream_length(1:number_of_realizations_per_distribution(G), 1:D-1, G)))) + sum(sum(sum(output_average_stream_length(1:rr, D, G)))) )/(nnz(output_average_stream_length(1:number_of_realizations_per_distribution(G), :, 1:G-1)) + nnz(output_average_stream_length(1:number_of_realizations_per_distribution(G), 1:D-1, G)) + nnz(output_average_stream_length(1:rr, D, G)) );
                cum_mean_average_stream_width(idx) = ( sum(sum(sum(output_average_stream_width(1:number_of_realizations_per_distribution(G), :, 1:G-1)))) + sum(sum(sum(output_average_stream_width(1:number_of_realizations_per_distribution(G), 1:D-1, G)))) + sum(sum(sum(output_average_stream_width(1:rr, D, G)))) )/(nnz(output_average_stream_width(1:number_of_realizations_per_distribution(G), :, 1:G-1)) + nnz(output_average_stream_width(1:number_of_realizations_per_distribution(G), 1:D-1, G)) + nnz(output_average_stream_width(1:rr, D, G)) );
                cum_mean_average_max_cell_x_location(idx) = ( sum(sum(sum(output_average_max_cell_x_location(1:number_of_realizations_per_distribution(G), :, 1:G-1)))) + sum(sum(sum(output_average_max_cell_x_location(1:number_of_realizations_per_distribution(G), 1:D-1, G)))) + sum(sum(sum(output_average_max_cell_x_location(1:rr, D, G)))) )/(nnz(output_average_max_cell_x_location(1:number_of_realizations_per_distribution(G), :,1:G-1)) + nnz(output_average_max_cell_x_location(1:number_of_realizations_per_distribution(G), 1:D-1, G)) + nnz(output_average_max_cell_x_location(1:rr, D, G)) );
                cum_mean_average_max_cell_y_width(idx) = ( sum(sum(sum(output_average_max_cell_y_width(1:number_of_realizations_per_distribution(G), :, 1:G-1)))) + sum(sum(sum(output_average_max_cell_y_width(1:number_of_realizations_per_distribution(G), 1:D-1, G)))) + sum(sum(sum(output_average_max_cell_y_width(1:rr, D, G)))) )/(nnz(output_average_max_cell_y_width(1:number_of_realizations_per_distribution(G), :,1:G-1)) + nnz(output_average_max_cell_y_width(1:number_of_realizations_per_distribution(G), 1:D-1, G)) + nnz(output_average_max_cell_y_width(1:rr, D, G)) );
                cum_mean_average_distance_to_nearest_cell(idx) = ( sum(sum(sum(output_average_distance_to_nearest_cell(1:number_of_realizations_per_distribution(G), :, 1:G-1)))) + sum(sum(sum(output_average_distance_to_nearest_cell(1:number_of_realizations_per_distribution(G), 1:D-1, G)))) + sum(sum(sum(output_average_distance_to_nearest_cell(1:rr, D, G)))) )/(nnz(output_average_distance_to_nearest_cell(1:number_of_realizations_per_distribution(G), :,1:G-1)) + nnz(output_average_distance_to_nearest_cell(1:number_of_realizations_per_distribution(G), 1:D-1, G)) + nnz(output_average_distance_to_nearest_cell(1:rr, D, G)) );
            end
            system(['rm -rf ',strcat(video_title,'_summary_statistics.mat')]);
            continue;
        end
    %% Compute all of the summary statistics you require:
        fprintf(['Group Number ',num2str(G), ' of ',num2str(number_of_groups),', Distribution Number ',num2str(D),' of ',num2str(number_of_distributions),'...\n']);
        Summary_Statistics = NCC_Model_Puncta_Multiple_Realizations_Consistency_Analysis(D, number_of_realizations_per_distribution(G));
        output_average_stream_length(1:number_of_realizations_per_distribution(G), D, G) = Summary_Statistics.Est_Stream_Length;
        output_average_stream_width(1:number_of_realizations_per_distribution(G), D, G) = Summary_Statistics.Est_Stream_Width;
        output_average_max_cell_x_location(1:number_of_realizations_per_distribution(G), D, G) = Summary_Statistics.Max_X_Length_of_Stream;
        output_average_max_cell_y_width(1:number_of_realizations_per_distribution(G), D, G) = Summary_Statistics.Max_Y_Length_of_Stream;
        output_average_distance_to_nearest_cell(1:number_of_realizations_per_distribution(G), D, G) = Summary_Statistics.Average_Min_Distance_per_Realization;
        
        for rr = 1:number_of_realizations_per_distribution(G)
            idx = number_of_distributions*sum(number_of_realizations_per_distribution(1:G-1)) + (D-1)*number_of_realizations_per_distribution(G) + rr;
            cum_mean_average_stream_length(idx) = ( sum(sum(sum(output_average_stream_length(1:number_of_realizations_per_distribution(G), :, 1:G-1)))) + sum(sum(sum(output_average_stream_length(1:number_of_realizations_per_distribution(G), 1:D-1, G)))) + sum(sum(sum(output_average_stream_length(1:rr, D, G)))) )/(nnz(output_average_stream_length(1:number_of_realizations_per_distribution(G), :, 1:G-1)) + nnz(output_average_stream_length(1:number_of_realizations_per_distribution(G), 1:D-1, G)) + nnz(output_average_stream_length(1:rr, D, G)) );
            cum_mean_average_stream_width(idx) = ( sum(sum(sum(output_average_stream_width(1:number_of_realizations_per_distribution(G), :, 1:G-1)))) + sum(sum(sum(output_average_stream_width(1:number_of_realizations_per_distribution(G), 1:D-1, G)))) + sum(sum(sum(output_average_stream_width(1:rr, D, G)))) )/(nnz(output_average_stream_width(1:number_of_realizations_per_distribution(G), :, 1:G-1)) + nnz(output_average_stream_width(1:number_of_realizations_per_distribution(G), 1:D-1, G)) + nnz(output_average_stream_width(1:rr, D, G)) );
            cum_mean_average_max_cell_x_location(idx) = ( sum(sum(sum(output_average_max_cell_x_location(1:number_of_realizations_per_distribution(G), :, 1:G-1)))) + sum(sum(sum(output_average_max_cell_x_location(1:number_of_realizations_per_distribution(G), 1:D-1, G)))) + sum(sum(sum(output_average_max_cell_x_location(1:rr, D, G)))) )/(nnz(output_average_max_cell_x_location(1:number_of_realizations_per_distribution(G), :,1:G-1)) + nnz(output_average_max_cell_x_location(1:number_of_realizations_per_distribution(G), 1:D-1, G)) + nnz(output_average_max_cell_x_location(1:rr, D, G)) );
            cum_mean_average_max_cell_y_width(idx) = ( sum(sum(sum(output_average_max_cell_y_width(1:number_of_realizations_per_distribution(G), :, 1:G-1)))) + sum(sum(sum(output_average_max_cell_y_width(1:number_of_realizations_per_distribution(G), 1:D-1, G)))) + sum(sum(sum(output_average_max_cell_y_width(1:rr, D, G)))) )/(nnz(output_average_max_cell_y_width(1:number_of_realizations_per_distribution(G), :,1:G-1)) + nnz(output_average_max_cell_y_width(1:number_of_realizations_per_distribution(G), 1:D-1, G)) + nnz(output_average_max_cell_y_width(1:rr, D, G)) );
            cum_mean_average_distance_to_nearest_cell(idx) = ( sum(sum(sum(output_average_distance_to_nearest_cell(1:number_of_realizations_per_distribution(G), :, 1:G-1)))) + sum(sum(sum(output_average_distance_to_nearest_cell(1:number_of_realizations_per_distribution(G), 1:D-1, G)))) + sum(sum(sum(output_average_distance_to_nearest_cell(1:rr, D, G)))) )/(nnz(output_average_distance_to_nearest_cell(1:number_of_realizations_per_distribution(G), :,1:G-1)) + nnz(output_average_distance_to_nearest_cell(1:number_of_realizations_per_distribution(G), 1:D-1, G)) + nnz(output_average_distance_to_nearest_cell(1:rr, D, G)) );
        end
    end
end

%% Compute A measure within each group and plot them:
A_measures_stream_length = zeros(number_of_distributions, number_of_groups);
A_measures_stream_width = zeros(number_of_distributions, number_of_groups);
A_measures_average_max_cell_x_location = zeros(number_of_distributions, number_of_groups);
A_measures_average_max_cell_y_width = zeros(number_of_distributions, number_of_groups);
A_measures_average_distance_to_nearest_cell = zeros(number_of_distributions, number_of_groups);

test_stream_length = zeros(number_of_distributions, number_of_groups, number_of_distributions);
test_stream_width = zeros(number_of_distributions, number_of_groups, number_of_distributions);
test_x = zeros(number_of_distributions, number_of_groups, number_of_distributions);
test_y = zeros(number_of_distributions, number_of_groups, number_of_distributions);
test_distance = zeros(number_of_distributions, number_of_groups, number_of_distributions);

for i = 1:number_of_groups
    for jj = 1:number_of_distributions
        for k = 1:number_of_distributions
            if jj == k
                continue;
            end
            test_stream_length(jj, i, k) =  getA_measure(output_average_stream_length(1:number_of_realizations_per_distribution(i), k,i), output_average_stream_length(1:number_of_realizations_per_distribution(i), jj ,i)) ;
            test_stream_width(jj, i, k) = getA_measure(output_average_stream_width(1:number_of_realizations_per_distribution(i), k,i), output_average_stream_width(1:number_of_realizations_per_distribution(i), jj ,i)) ;
            test_x(jj, i, k) = getA_measure(output_average_max_cell_x_location(1:number_of_realizations_per_distribution(i), k,i), output_average_max_cell_x_location(1:number_of_realizations_per_distribution(i), jj ,i)) ;
            test_y(jj, i, k) = getA_measure(output_average_max_cell_y_width(1:number_of_realizations_per_distribution(i), k,i), output_average_max_cell_y_width(1:number_of_realizations_per_distribution(i), jj ,i)) ;
            test_distance(jj, i, k) = getA_measure(output_average_distance_to_nearest_cell(1:number_of_realizations_per_distribution(i), k,i), output_average_distance_to_nearest_cell(1:number_of_realizations_per_distribution(i), jj ,i)) ;
        end
    end
end
for i = 1:number_of_groups
    for jj = 1:number_of_distributions
        if jj == 1
            continue;
        end
        A_measures_stream_length(jj, i) = getA_measure(output_average_stream_length(1:number_of_realizations_per_distribution(i), 1,i), output_average_stream_length(1:number_of_realizations_per_distribution(i), jj ,i)) ;
        A_measures_stream_width(jj, i) = getA_measure(output_average_stream_width(1:number_of_realizations_per_distribution(i), 1,i), output_average_stream_width(1:number_of_realizations_per_distribution(i), jj ,i)) ;
        A_measures_average_max_cell_x_location(jj, i) = getA_measure(output_average_max_cell_x_location(1:number_of_realizations_per_distribution(i), 1,i), output_average_max_cell_x_location(1:number_of_realizations_per_distribution(i), jj ,i)) ;
        A_measures_average_max_cell_y_width(jj, i) = getA_measure(output_average_max_cell_y_width(1:number_of_realizations_per_distribution(i), 1,i), output_average_max_cell_y_width(1:number_of_realizations_per_distribution(i), jj ,i)) ;
        A_measures_average_distance_to_nearest_cell(jj, i) = getA_measure(output_average_distance_to_nearest_cell(1:number_of_realizations_per_distribution(i), 1,i), output_average_distance_to_nearest_cell(1:number_of_realizations_per_distribution(i), jj ,i)) ;
    end % for jj
    fig = figure('units', 'inches', 'position', [0,0,7,5]);
    p1 = plot(2:number_of_distributions, A_measures_stream_length(2:number_of_distributions,i), '-*');
    hold on;
    p2 = plot(2:number_of_distributions, A_measures_stream_width(2:number_of_distributions,i), '-*');
    p3 = plot(2:number_of_distributions, A_measures_average_max_cell_x_location(2:number_of_distributions,i), '-*');
    p4 = plot(2:number_of_distributions, A_measures_average_max_cell_y_width(2:number_of_distributions,i), '-*');
    p5 = plot(2:number_of_distributions, A_measures_average_distance_to_nearest_cell(2:number_of_distributions,i), '-*');
    xlabel('Group Number');
    ylabel('A-measures');
    legend([p1(1), p2(1), p3(1), p4(1), p5(1)], 'Stream Length', 'Stream Width', 'Average x-location', 'Average y-width', 'Average nearest distance');
    saveas(fig, strcat(save_path,'A_measures_Group_', num2str(i),'.png'), 'png');
    saveas(fig, strcat(save_path,'A_measures_Group_', num2str(i),'.svg'), 'svg');
end % for i
    fig2 = figure('units', 'inches', 'position', [0,0,7,5]);
    subplot(5,1,1);
    plot(1:length(cum_mean_average_stream_length), cum_mean_average_stream_length);
    ylabel('Stream Length, Cum. Mean');
    subplot(5,1,2);
    plot(1:length(cum_mean_average_stream_width), cum_mean_average_stream_width);
    ylabel('Stream Width, Cum. Mean')
    subplot(5,1,3);
    plot(1:length(cum_mean_average_max_cell_x_location), cum_mean_average_max_cell_x_location);
    ylabel('Max X Location, Cum. Mean');
    subplot(5,1,4);
    plot(1:length(cum_mean_average_max_cell_y_width), cum_mean_average_max_cell_y_width);
    ylabel('Max Y Width, Cum. Mean');
    subplot(5,1,5);
    plot(1:length(cum_mean_average_distance_to_nearest_cell), cum_mean_average_distance_to_nearest_cell);
    xlabel('Total Number of Realizations Run')
    ylabel('Average Nearest Distance, Cum. Mean')
    saveas(fig2, strcat(save_path,'Cum_Mean_Plots.png'), 'png');
    saveas(fig2, strcat(save_path,'Cum_Mean_Plots.svg'), 'svg');
    
[max_A_measure_stream_length, max_A_index_stream_length] = max(A_measures_stream_length, [], 1);
[max_A_measure_stream_width, max_A_index_stream_width] = max(A_measures_stream_width, [], 1);
[max_A_measure_average_max_cell_x_location, max_A_index_average_max_cell_x_location] = max(A_measures_average_max_cell_x_location, [], 1);
[max_A_measure_average_max_cell_y_width, max_A_index_average_max_cell_y_width] = max(A_measures_average_max_cell_y_width, [], 1);
[max_A_measure_average_distance_to_nearest_cell, max_A_index_average_distance_to_nearest_cell] = max(A_measures_average_distance_to_nearest_cell, [], 1);

test_alt_A_measure_stream_length = max( max(test_stream_length,[], 1), [], 3);
test_alt_A_measure_stream_width = max( max(test_stream_width,[], 1), [], 3);
test_alt_A_measure_average_max_cell_x_location = max( max(test_x,[], 1), [], 3);
test_alt_A_measure_average_max_cell_y_width = max( max(test_y,[], 1), [], 3);
test_alt_A_measure_average_distance_to_nearest_cell = max( max(test_distance,[], 1), [], 3);

threshold_small = 0.56;
threshold_medium = 0.64;
threshold_large = 0.71;

best_index_stream_length = find(max_A_measure_stream_length < threshold_small, 1, 'first');
best_index_stream_width = find(max_A_measure_stream_width < threshold_small, 1, 'first');
best_index_average_max_cell_x_location = find(max_A_measure_average_max_cell_x_location < threshold_small, 1, 'first');
best_index_average_max_cell_y_width = find(max_A_measure_average_max_cell_y_width < threshold_small, 1, 'first');
best_index_average_distance_to_nearest_cell = find(max_A_measure_average_distance_to_nearest_cell < threshold_small, 1, 'first');

test_best_index_stream_length = find(test_alt_A_measure_stream_length < threshold_small, 1, 'first');
test_best_index_stream_width = find(test_alt_A_measure_stream_width < threshold_small, 1, 'first');
test_best_index_average_max_cell_x_location = find(test_alt_A_measure_average_max_cell_x_location < threshold_small, 1, 'first');
test_best_index_average_max_cell_y_width = find(test_alt_A_measure_average_max_cell_y_width < threshold_small, 1, 'first');
test_best_index_average_distance_to_nearest_cell = find(test_alt_A_measure_average_distance_to_nearest_cell < threshold_small, 1, 'first');

Optimal_Realizations = {{'best_number_stream_length', 'best_number_stream_width', 'best_number_average_max_cell_x_location', 'best_number_average_max_cell_y_width', 'best_number_average_distance_to_nearest_cell'}; {number_of_realizations_per_distribution(best_index_stream_length), number_of_realizations_per_distribution(best_index_stream_width), number_of_realizations_per_distribution(best_index_average_max_cell_x_location), number_of_realizations_per_distribution(best_index_average_max_cell_y_width), number_of_realizations_per_distribution(best_index_average_distance_to_nearest_cell)}};
Test_Optimal_Realizations = {{'best_number_stream_length', 'best_number_stream_width', 'best_number_average_max_cell_x_location', 'best_number_average_max_cell_y_width', 'best_number_average_distance_to_nearest_cell'};{number_of_realizations_per_distribution(test_best_index_stream_length), number_of_realizations_per_distribution(test_best_index_stream_width), number_of_realizations_per_distribution(test_best_index_average_max_cell_x_location), number_of_realizations_per_distribution(test_best_index_average_max_cell_y_width), number_of_realizations_per_distribution(test_best_index_average_distance_to_nearest_cell)}};
save(strcat(save_path, 'optimal_number_of_realizations.mat'), 'Optimal_Realizations', 'Test_Optimal_Realizations');

fig = figure('units', 'inches', 'position', [0,0,7,5]);
p1 = plot(number_of_realizations_per_distribution, max_A_measure_stream_length, '-*');
hold on;
p2 = plot(number_of_realizations_per_distribution, max_A_measure_stream_width, '-*');
p3 = plot(number_of_realizations_per_distribution, max_A_measure_average_max_cell_x_location, '-*');
p4 = plot(number_of_realizations_per_distribution, max_A_measure_average_max_cell_y_width, '-*');
p5 = plot(number_of_realizations_per_distribution, max_A_measure_average_distance_to_nearest_cell, '-*');
xlabel('Number of Realizations');
ylabel('A_measure');
legend([p1(1), p2(1), p3(1), p4(1), p5(1)], 'Stream Length', 'Stream Width', 'Average x-location', 'Average y-width', 'Average nearest distance');
fplot(@(x) 0*x + threshold_small, [number_of_realizations_per_distribution(1),number_of_realizations_per_distribution(end)] , '--g');
fplot(@(x) 0*x + threshold_medium, [number_of_realizations_per_distribution(1),number_of_realizations_per_distribution(end)], '--y');
fplot(@(x) 0*x + threshold_large, [number_of_realizations_per_distribution(1),number_of_realizations_per_distribution(end)], '--r');
saveas(fig, strcat(save_path,'Max_A_measures.png'), 'png');
saveas(fig, strcat(save_path,'Max_A_measures.svg'), 'svg');

fig = figure('units', 'inches', 'position', [0,0,7,5]);
p1 = plot(number_of_realizations_per_distribution, test_alt_A_measure_stream_length, '-*');
hold on;
p2 = plot(number_of_realizations_per_distribution, test_alt_A_measure_stream_width, '-*');
p3 = plot(number_of_realizations_per_distribution, test_alt_A_measure_average_max_cell_x_location, '-*');
p4 = plot(number_of_realizations_per_distribution, test_alt_A_measure_average_max_cell_y_width, '-*');
p5 = plot(number_of_realizations_per_distribution, test_alt_A_measure_average_distance_to_nearest_cell, '-*');
xlabel('Number of Realizations');
ylabel('A_measure');
legend([p1(1), p2(1), p3(1), p4(1), p5(1)], 'Stream Length', 'Stream Width', 'Average x-location', 'Average y-width', 'Average nearest distance');
fplot(@(x) 0*x + threshold_small, [number_of_realizations_per_distribution(1),number_of_realizations_per_distribution(end)] , '--g');
fplot(@(x) 0*x + threshold_medium, [number_of_realizations_per_distribution(1),number_of_realizations_per_distribution(end)], '--y');
fplot(@(x) 0*x + threshold_large, [number_of_realizations_per_distribution(1),number_of_realizations_per_distribution(end)], '--r');
saveas(fig, strcat(save_path,'Max_Test_A_measures.png'), 'png');
saveas(fig, strcat(save_path,'Max_Test_A_measures.svg'), 'svg');

%% Supplementary Functions
function [A_measure, scaled_A_measure] = getA_measure(x0, x1) 
[~,~,stats] = ranksum(x0,x1); 
% Compute the A measure 
A_measure=(stats.ranksum/length(x0) - (length(x0)+1)/2)/length(x1); 
% Compute the scaled A measure 
scaled_A_measure=0.5+abs(0.5 - A_measure);
end

% function [A_measure, scaled_A_measure] = getA_measure_naive(x0, x1)
% % Compute the A measure 
% A_measure = 0; 
% for i = 1:length(x0) 
%     for j = 1:length(x1) 
%         if(x0(i)>x1(j)) 
%             A_measure = A_measure + 1;
%         elseif(x0(i)==x1(j)) 
%             A_measure = A_measure + 0.5;
%         elseif(x0(i)<x1(j))
%             A_measure = A_measure + 0;
%         end
%     end
% end
% A_measure = A_measure/(length(x0)*length(x1)); 
% % Compute the scaled A measure 
% if(A_measure>=0.5) 
%     scaled_A_measure = A_measure;
% else
%     scaled_A_measure = 1-A_measure;
% end
% end