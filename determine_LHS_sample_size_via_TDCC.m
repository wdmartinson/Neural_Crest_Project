%% Determine Optimal LHS Sample Size, via TDCCs, using the procedure outlined in Marino et al. (2008)
% Note that we are subsampling from the 1000 parameter regimes already
% computed, to save ourselves extra computational time.
compute_lhs_sample_size = true;
if compute_lhs_sample_size
rng(3); % Set the random seed, for reproducability
N = cell(6, 1);
num_replicates = 30;
sample_sizes = [100 200 300 400 500 1000];
for ii = 1:6
    N{ii} = zeros(sample_sizes(ii), num_replicates);
    for jj = 1:num_replicates
        N{ii}(:,jj) = randperm(1000, sample_sizes(ii)); % Average the PRCCs over 100 subsamples
    end
end
% N{1} = randi(1000, 100, num_replicates); % Average the PRCCs over 100 subsamples
% N{2} = randi(1000, 200, num_replicates);
% N{3} = randi(1000, 300, num_replicates);
% N{4} = randi(1000, 400, num_replicates);
% N{5} = randi(1000, 500, num_replicates);
% N{6} = repmat((1:1000)', 1, num_replicates);
num_sample_sizes = length(N); % number of sample sizes (N = 100, N = 200, N = 300, N = 400, N = 500, N = 1000)

% Load all of the output statistics
load('/Volumes/easystore/NCC_puncta_model_output_data/PRCC_data.mat', 'output_average_max_cell_x_location', 'output_average_max_cell_y_width', 'output_average_distance_to_nearest_cell');
load('/Volumes/easystore/NCC_puncta_model_output_data/LHS_output_extra_stats.mat', 'output_average_follower_order_parameter');
num_outputs = 4;
% Load the parameter regime matrix:
LHS_Parameter_Matrix = readmatrix('/Volumes/easystore/NCC_puncta_model_output_data/LHS_Parameter_Values_Matrix.csv');
num_parameters = size(LHS_Parameter_Matrix, 2); % number of parameters
Names = { 'fibronectin_puncta_wavelength',...
                     'filopodia_sensing_radius',...
                     'bias_towards_VM_direction',...
                     'average_time_until_next_filopodia_drop',...
                     'half_life_of_puncta_orientation_change',...
                     'y_length_of_cell_entrance_strip',...
                     'default_cell_radius',...
                     'default_cell_speed',...
                     'followers_cell_cell_repulsion_strength'};

PRCC_Matrix = zeros(num_sample_sizes, num_parameters, num_outputs);
temp = zeros(num_replicates, num_parameters, 4, num_sample_sizes);
for ii = 1:num_sample_sizes
    for pp = 1:num_replicates
        idx = N{ii}(:, pp);
        regimes = LHS_Parameter_Matrix(idx, :);
        temp(pp, :, 1, ii) = PRCC(regimes, output_average_max_cell_x_location(idx), 720, Names, 0.01); % PRCC value for this
        temp(pp, :, 2, ii) = PRCC(regimes, output_average_max_cell_y_width(idx), 720, Names, 0.01); % PRCC value for this
        temp(pp, :, 3, ii) = PRCC(regimes, output_average_distance_to_nearest_cell(idx), 720, Names, 0.01); % PRCC value for this
        temp(pp, :, 4, ii) = PRCC(regimes, output_average_follower_order_parameter(idx), 720, Names, 0.01); % PRCC value for this
    end % for pp
%     PRCC_Matrix(ii, :, 1) = PRCC(regimes, output_average_max_cell_x_location(idx), 720, Names, 0.01); % PRCC value for this
%     PRCC_Matrix(ii, :, 2) = PRCC(regimes, output_average_max_cell_y_width(idx), 720, Names, 0.01); % PRCC value for this
%     PRCC_Matrix(ii, :, 3) = PRCC(regimes, output_average_distance_to_nearest_cell(idx), 720, Names, 0.01); % PRCC value for this
%     PRCC_Matrix(ii, :, 4) = PRCC(regimes, output_average_follower_order_parameter(idx), 720, Names, 0.01); % PRCC value for this
    for zz = 1:num_outputs
        PRCC_Matrix(ii, :, zz) = mean(temp(:,:, zz, ii), 1);
    end % for zz
end % for ii

TDCCs = zeros(num_sample_sizes-1, num_outputs);
significant_TDCCs = TDCCs;

for rr = 1:num_outputs
%     [Tnew,~]=ranking(PRCC_Matrix(:,:,rr)');
    [Tnew,~]=ranking(abs(PRCC_Matrix(:,:,rr))');
    for A = 1:num_sample_sizes-1
        [~, TDCCs(A, rr), significant_TDCCs(A, rr)]=topdowncorr(Tnew(:, A:A+1));
    end
end
significant_TDCCs_logical = (significant_TDCCs < 0.01); % Determine which TDCCs have significant p-values


else % compute the eFAST sample size instead
%% Optimal eFAST sample size, using TDCCs and the method of Marino et al. (2008)
% Note that we are taking subsamples rather than computing new datasets
% Also, we only conduct our analysis on the TDCCs (as it would otherwise be
% impossible to compute Sobol Indices from subsamples at these frequencies)
rng(4); % Reset random seed, for reproducability
N = cell(5, 1);
num_replicates = 100;
N{1} = randi(131*4*3, 65*4, num_replicates); % Average the PRCCs over 100 subsamples
N{2} = randi(131*4*3, 105*4, num_replicates);
N{3} = randi(131*4*3, 65*4*2, num_replicates);
N{4} = randi(131*4*3, 105*4*2, num_replicates);
N{5} = randi(131*4*3, 65*4*3, num_replicates);
N{6} = randi(131*4*3, 105*4*3, num_replicates);
N{7} = repmat((1:131*4*3)', 1, num_replicates);
num_sample_sizes = length(N); % number of sample sizes (N = 100, N = 200, N = 300, N = 400, N = 500, N = 1000)

load('/Volumes/easystore/NCC_puncta_model_output_data_eFAST/eFAST_output_stats.mat', 'output_average_max_cell_x_location', 'output_average_max_cell_y_width', 'output_average_distance_to_nearest_cell');
load('/Volumes/easystore/NCC_puncta_model_output_data_eFAST/eFAST_output_extra_stats.mat', 'output_average_follower_order_parameter');
load('/Volumes/easystore/NCC_puncta_model_output_data_eFAST/eFAST_sample_matrix.mat', 'eFAST_sample_matrix');
num_outputs = 4;
num_parameters = size(eFAST_sample_matrix, 2);
Names = {'filopodia_sensing_radius',...
                     'bias_towards_VM_direction',...
                     'followers_cell_cell_repulsion_strength',...
                     'dummy_parameter'};
output_average_max_cell_x_location = reshape(output_average_max_cell_x_location, size(eFAST_sample_matrix,1)*size(eFAST_sample_matrix,3), 1 )';
output_average_max_cell_y_width = reshape(output_average_max_cell_y_width, size(eFAST_sample_matrix,1)*size(eFAST_sample_matrix,3), 1 )';
output_average_distance_to_nearest_cell = reshape(output_average_distance_to_nearest_cell, size(eFAST_sample_matrix,1)*size(eFAST_sample_matrix,3), 1 )';
output_average_follower_order_parameter = reshape(output_average_follower_order_parameter, size(eFAST_sample_matrix,1)*size(eFAST_sample_matrix,3), 1 )';
eFAST_sample_matrix = [eFAST_sample_matrix(:, :, 1); eFAST_sample_matrix(:, :, 2); eFAST_sample_matrix(:, :, 3)];
                 
PRCC_Matrix = zeros(num_sample_sizes, num_parameters, num_outputs);
temp = zeros(num_replicates, num_parameters, 4, num_sample_sizes);
for ii = 1:num_sample_sizes
    for pp = 1:num_replicates
        idx = N{ii}(:, pp); % Take all of these together fo
        % Concatenate the replicate samples together?
        regimes = eFAST_sample_matrix(idx, :);
        temp(pp, :, 1, ii) = PRCC(regimes, output_average_max_cell_x_location(idx), 720, Names, 0.01); % PRCC value for this
        temp(pp, :, 2, ii) = PRCC(regimes, output_average_max_cell_y_width(idx), 720, Names, 0.01); % PRCC value for this
        temp(pp, :, 3, ii) = PRCC(regimes, output_average_distance_to_nearest_cell(idx), 720, Names, 0.01); % PRCC value for this
        temp(pp, :, 4, ii) = PRCC(regimes, output_average_follower_order_parameter(idx), 720, Names, 0.01); % PRCC value for this
    end % for pp
    for zz = 1:num_outputs
        PRCC_Matrix(ii, :, zz) = mean(temp(:,:, zz, ii), 1);
    end % for zz
end % for ii                

TDCCs = zeros(num_sample_sizes-1, num_outputs);
significant_TDCCs = TDCCs;
for rr = 1:num_outputs
    [Tnew, ~] = ranking(abs(PRCC_Matrix(:,:,rr))');
    for A = 1:num_sample_sizes-1
        [~, TDCCs(A, rr), significant_TDCCs(A, rr)]=topdowncorr(Tnew(:, A:A+1));
    end
end
significant_TDCCs_logical = (significant_TDCCs < 0.01); % Determine which TDCCs have significant p-values

end % if compute_lhs_sample_size
