% main_puncta_experiments.m : Code that performs in silico versions of
% experiments suggested by the Kulesa lab.
%% Make sure you are using the correct version of Python and g++
str = computer;
desktop = true;
save_path = './NCC_puncta_model_output_data_experiments/';
if str(1) == 'M'
    desktop = false;
    setenv('PATH', '/usr/local/bin:/usr/local/sbin:/usr/bin:/bin:/usr/sbin:/sbin');
    save_path = './NCC_puncta_model_output_data_experiments/';
end
%% Create cell with all of the parameter values you wish to change
parameters.names =  {'scenario',...
                     'fibronectin_puncta_wavelength',...
                     'bias_towards_VM_direction',...
                     'average_time_until_next_filopodia_drop',...
                     'y_length_of_cell_entrance_strip',...
                     'followers_cell_cell_repulsion_strength'};
                 
%% Create matrix with all of the parameter values 
% Note that not all of the parameter values will change with each
% experiment.
number_of_parameters = length(parameters.names);

% Write the parameter value matrix as a CSV file for safe keeping:
if ~isfile(strcat(save_path, 'experiments_parameter_matrix.mat'))
    % 11 Scenarios to run, excluding rescue experiments
    parameters.values = zeros(12, number_of_parameters);
    parameters.values(:,1) = (1:12)'; % scenario numbers
    parameters.values(:,2) = 20*ones(12,1); % fibronectin_puncta_wavelength
    parameters.values(:,3) = 0.5*ones(12,1); % bias_towards_VM_direction
    parameters.values(:,4) = 30*ones(12, 1); % average_time_until_next_filopodia_drop
    parameters.values(:,5) = 120*ones(12, 1); % y_length_of_cell_entrance_strip
    parameters.values(:,6) = 0.5*ones(12, 1); % cell_cell_repulsion_strength (applies to leaders and followers)
    
    % Scenario B: Test what happens with no cell-cell repulsion:
    parameters.values(2, end) = 0;
    % Scenario C: Test what happens when the FN puncta spacing is very
    %             large (over 2x the filopodia radius of the cell)
    parameters.values(3, 2) = 60; % microns
    % Scenario J: Significantly reduce the number of leaders by decreasing
    % the size of the cell entrance strip:
    parameters.values(10, 5) = 30; % microns
    
    % Rescue experiments:
    % For the case in which c_i = 0:
    rescue_2_1 = parameters.values(2,:);
    rescue_2_1(2) = 15; % fibronectin_puncta_wavelength (um)
    rescue_2_2 = parameters.values(2,:);
    rescue_2_2(3) = 0.25; % bias_towards_VM_direction
    rescue_2_3 = parameters.values(2,:);
    rescue_2_3(4) = 10; % average_time_until_next_filopodia_drop (min)
    rescue_2_4 = parameters.values(2,:);
    rescue_2_4(end) = 1e-5; % cell_cell_repulsion_strength (applies to leaders and followers)
    
    overexpression_2_5 = parameters.values(2,:);
    overexpression_2_5(end) = 2; % cell_cell_repulsion_strength (applies to leaders and followers)
    
    % For the case in which R_filo - lambda_FN << 0
    rescue_3_1 = parameters.values(3,:);
    rescue_3_1(3) = 0.25; % bias_towards_VM_direction
    rescue_3_2 = parameters.values(3,:);
    rescue_3_2(4) = 10; % average_time_until_next_filopodia_drop (min)
    rescue_3_3 = parameters.values(3,:);
    rescue_3_3(end) = 2; % cell_cell_repulsion_strength (applies to leaders and followers)
    
    parameters.values = [parameters.values; rescue_2_1; rescue_2_2; ...
        rescue_2_3; rescue_2_4; overexpression_2_5; rescue_3_1; rescue_3_2; rescue_3_3];
    
    % Latin Hypercube sampling, to see how results are affected by
    % perturbations in metaparameters (lambda_FN, rho, and c_i):
    rng(1); % for reproducability
    lhs_mat = lhsdesign(30, 3);
    control_perturbations = repmat(parameters.values(1,:), 30, 1);
    control_perturbations(:, 2) = 10 + 20*lhs_mat(:,1); % fibronectin_puncta_wavelength in [10, 30]
    control_perturbations(:, 3) = lhs_mat(:,2); % bias_toward_VM_direction in [0, 1]
    control_perturbations(:, end) = 2*lhs_mat(:, end); % cell_cell_repulsion_strength (applies to leaders and followers)
    
    parameters.values = [parameters.values; control_perturbations];
    
    % Finally, run similar experiments to Scenario 2 and Scenario 8-9, 12, 
    % in which the contact guidance cue is much less pronounced.
    less_contact_guidance_expts = repmat(parameters.values(2,:), 5, 1);
    less_contact_guidance_expts(:,3) = 0.95*ones(size(less_contact_guidance_expts, 1), 1); % bias_toward_VM_direction
    less_contact_guidance_expts(2, 2) = 15; % fibronectin_puncta_spacing
    less_contact_guidance_expts(3, 4) = 10; % average_time_until_next_filopodia_drop (min)
    less_contact_guidance_expts(4, end) = 1e-5; % cell_cell_repulsion_strength (applies to leaders and followers)
    less_contact_guidance_expts(5, end) = 2; % cell_cell_repulsion_strength (applies to leaders and followers)
    
    less_contact_guidance_expts = [less_contact_guidance_expts; parameters.values(8:9, :); parameters.values(12,:)];
    less_contact_guidance_expts(:,3) = 0.95*ones(size(less_contact_guidance_expts, 1), 1); % bias_toward_VM_direction
    
    parameters.values = [parameters.values; less_contact_guidance_expts];
    
    plateau_effect_expts = parameters.values(1,:); 
    % Control expt for 
    plateau_effect_expts(1, 3) = 0.95; % bias_toward_VM_direction
    plateau_effect_expts = [plateau_effect_expts; repmat(parameters.values(13,:), 2, 1); repmat(parameters.values(14,:), 2, 1); repmat(parameters.values(15,:), 2, 1); parameters.values(17,:); repmat(parameters.values(18,:), 2, 1); repmat(parameters.values(19,:), 2, 1); parameters.values(20,:)];
    plateau_effect_expts(2, 2) = 10; % fibronectin_puncta_wavelength
    plateau_effect_expts(3, 2) = 5; % fibronectin_puncta_wavelength
    plateau_effect_expts(4, 3) = 0.1; % bias_toward_VM_direction
    plateau_effect_expts(5, 3) = 0.05; % bias_toward_VM_direction
    plateau_effect_expts(6, 4) = 5; % average_time_until_next_filopodia_drop (min)
    plateau_effect_expts(7, 4) = 1; % average_time_until_next_filopodia_drop (min)
    plateau_effect_expts(8, end) = 5; % cell_cell_repulsion_strength (applies to leaders and followers)
    plateau_effect_expts(9, 3) = 0.1; % bias_toward_VM_direction
    plateau_effect_expts(10, 3) = 0.05; % bias_toward_VM_direction
    plateau_effect_expts(11, 4) = 5; % average_time_until_next_filopodia_drop (min)
    plateau_effect_expts(12, 4) = 1; % average_time_until_next_filopodia_drop (min)
    plateau_effect_expts(13, end) = 5; % cell_cell_repulsion_strength (applies to leaders and followers)
    
    parameters.values = [parameters.values; plateau_effect_expts];
    
    new_experiments = repmat(parameters.values(1,:), 12, 1);
    new_experiments(1,1) = 13; % scenario
    new_experiments(2,1) = 14; % scenario
    new_experiments(3,1) = 15; % scenario
    new_experiments(4:7, 1) = 16*ones(4,1); % scenario
    new_experiments(4, 4) = 90; % average_time_until_next_filopodia_drop (min)
    new_experiments(5, 4) = 30; % average_time_until_next_filopodia_drop (min)
    new_experiments(6, 4) = 10; % average_time_until_next_filopodia_drop (min)
    new_experiments(7, 4) = 5; % average_time_until_next_filopodia_drop (min)
    new_experiments(8, 4) = 10; % average_time_until_next_filopodia_drop (min)
    new_experiments(9, 1) = 17; % scenario
    new_experiments(9, 4) = 10; % average_time_until_next_filopodia_drop (min)
    new_experiments(10, 1) = 18; % scenario
    new_experiments(11, 1) = 19; % scenario
    new_experiments(12, 1) = 17; % scenario
    
    parameters.values = [parameters.values; new_experiments];
    
    secretionless_expts = [parameters.values(1,:); parameters.values(59, :); control_perturbations; repmat(parameters.values(1,:), 4, 1)];
    secretionless_expts(:, 1) = ones(size(secretionless_expts, 1), 1).*20; % scenario
    secretionless_expts(end-3, end) = 2; % cell_cell_repulsion_strength (applies to leaders and followers)
    secretionless_expts(end-2, end) = 5; % cell_cell_repulsion_strength (applies to leaders and followers)
    secretionless_expts(end-1, 2) = 10; % fibronectin_puncta_wavelength
    secretionless_expts(end, 2) = 60; % fibronectin_puncta_wavelength
    
    parameters.values = [parameters.values; secretionless_expts];
    
    overexpression_expts = repmat(parameters.values(1,:), 2, 1);
    overexpression_expts(:, 1) = ones(2,1).*21; % scenario
    overexpression_expts(1, 4) = 10; % average_time_until_next_filopodia_drop (min)
    overexpression_expts(2, 4) = 5; % average_time_until_next_filopodia_drop (min)
    
    parameters.values = [parameters.values; overexpression_expts];
    
    strip_ICs = repmat(parameters.values(1,:), 14, 1);
    strip_ICs(1, 1) = 22; % scenario
    strip_ICs(2, 1) = 23; % scenario
    strip_ICs(3, 1) = 24; % scenario
    strip_ICs(4, 1) = 25; % scenario
    strip_ICs(5, 1) = 26; % scenario
    strip_ICs(6:11, 1) = ones(length(6:11), 1).*23; % scenario
    strip_ICs(6, end) = 2; % cell_cell_repulsion_strength (applies to leaders and followers)
    strip_ICs(7, end) = 5; % cell_cell_repulsion_strength (applies to leaders and followers)
    strip_ICs(8, 2) = 10; % fibronectin_puncta_wavelength
    strip_ICs(9, 4) = 10; % average_time_until_next_filopodia_drop (min)
    strip_ICs(10, 2) = 60; % fibronectin_puncta_wavelength
    strip_ICs(11, 2) = 60; % fibronectin_puncta_wavelength
    strip_ICs(11, 4) = 10; % average_time_until_next_filopodia_drop (min)
    strip_ICs(12, 1) = 27; % scenario
    strip_ICs(13, 1) = 28; % scenario
    strip_ICs(14, 1) = 28; % scenario
    strip_ICs(14, end) = 2; % cell_cell_repulsion_strength (applies to leaders and followers)
    
    parameters.values = [parameters.values; strip_ICs];
    
    chemotaxis_regimes = repmat(parameters.values(1,:), 24, 1);
    chemotaxis_regimes(7:12, 1) = 15*ones(6, 1); % scenario
    chemotaxis_regimes(13:18, 1) = 20*ones(6, 1); % scenario
    chemotaxis_regimes(19:24, 2) = 60*ones(6, 1); % fibronectin_puncta_wavelength
    chemotaxis_regimes(22:24, 4) = 10*ones(3, 1); % average_time_until_next_filopodia_drop (min)
    
    Chemotaxis_Parameters = repmat([1, 0, 0.25], 3, 1);
    Chemotaxis_Parameters(2, 3) = 0.5; Chemotaxis_Parameters(3, 3) = 0.75;
    Chemotaxis_Parameters = repmat(Chemotaxis_Parameters, 2, 1);
    Chemotaxis_Parameters(4:end, 2) = ones(3, 1);
    Chemotaxis_Parameters = repmat(Chemotaxis_Parameters, 3, 1);
    Chemotaxis_Parameters = [Chemotaxis_Parameters; repmat(Chemotaxis_Parameters(1:3, :), 2, 1)];
    
    parameters.values = [parameters.values; chemotaxis_regimes];
    
    Stream_Break_Expts = repmat(parameters.values(1,:), 22, 1);
    Stream_Break_Expts(3, 3) = 0.33; % bias_toward_VM_direction
    Stream_Break_Expts(4, 4) = 10; % average_time_until_next_filopodia_drop (min)
    Stream_Break_Expts(5, end) = 0.05; % cell_cell_repulsion_strength (applies to leaders and followers)
    Stream_Break_Expts(6, end) = 2; % cell_cell_repulsion_strength (applies to leaders and followers)
    Stream_Break_Expts(7:9, 1) = (29:31)'; % scenario
    Stream_Break_Expts(10:15, 3) = 0.33*ones(6,1); % bias_toward_VM_direction
    Stream_Break_Expts(10, 4) = 10; % average_time_until_next_filopodia_drop (min)
    Stream_Break_Expts(11, 1) = 29; % scenario
    Stream_Break_Expts(12, end) = 0.05; % cell_cell_repulsion_strength (applies to leaders and followers)
    Stream_Break_Expts(13, end) = 2; % cell_cell_repulsion_strength (applies to leaders and followers)
    Stream_Break_Expts(14, 1) = 30; % scenario
    Stream_Break_Expts(15, 1) = 31; % scenario
    Stream_Break_Expts(16, 1) = 29; % scenario
    Stream_Break_Expts(16:18, 4) = 10*ones(3,1); % average_time_until_next_filopodia_drop (min)
    Stream_Break_Expts(17, end) = 0.05; % cell_cell_repulsion_strength (applies to leaders and followers)
    Stream_Break_Expts(18, end) = 2; % cell_cell_repulsion_strength (applies to leaders and followers)
    Stream_Break_Expts(19, end) = 0.05; % cell_cell_repulsion_strength (applies to leaders and followers)
    Stream_Break_Expts(19, 1) = 30; % scenario
    Stream_Break_Expts(20, end) = 2; % cell_cell_repulsion_strength (applies to leaders and followers)
    Stream_Break_Expts(20, 1) = 32; % scenario
    Stream_Break_Expts(21, 4) = 10; % average_time_until_next_filopodia_drop (min)
    Stream_Break_Expts(21, end) = 0.05; % cell_cell_repulsion_strength (applies to leaders and followers)
    Stream_Break_Expts(21, 1) = 30; % scenario
    Stream_Break_Expts(22, 3) = 0.33; % bias_toward_VM_direction
    Stream_Break_Expts(22, 4) = 10; % average_time_until_next_filopodia_drop (min)
    Stream_Break_Expts(22, end) = 0.05; % cell_cell_repulsion_strength (applies to leaders and followers)

    Stream_Break_Expts = repmat(Stream_Break_Expts, 3, 1);
    Stream_Break_Expts(23:44, 2) = 60*ones(22, 1); % fibronectin_puncta_wavelength
    
    Stream_Break_Expts(45:end, 1) = 20*ones(22, 1); % scenario
    Stream_Break_Expts(48, 1) = 33; % scenario
    Stream_Break_Expts(51, 1) = 34; % scenario
    Stream_Break_Expts(52, 1) = 35; % scenario
    Stream_Break_Expts(53, 1) = 36; % scenario
    Stream_Break_Expts(54, 1) = 33; % scenario
    Stream_Break_Expts(55, 1) = 34; % scenario
    Stream_Break_Expts(58, 1) = 35; % scenario
    Stream_Break_Expts(59, 1) = 36; % scenario
    Stream_Break_Expts(60, 1) = 37; % scenario
    Stream_Break_Expts(61:62, 1) = 33*ones(2,1); % scenario
    Stream_Break_Expts(63, 1) = 35; % scenario
    Stream_Break_Expts(64, 1) = 38; % scenario
    Stream_Break_Expts(65, 1) = 39; % scenario
    Stream_Break_Expts(66, 1) = 33; % scenario
    
    parameters.values = [parameters.values; Stream_Break_Expts];
    
    Other_Rescue_Expts = repmat(parameters.values(1,:), 20, 1);
    Other_Rescue_Expts(1:5,1) = 40*ones(5,1); % scenario
    Other_Rescue_Expts(2, 4) = 10; % average_time_until_next_filopodia_drop (min)
    Other_Rescue_Expts(3, 1) = 41; %scenario
    Other_Rescue_Expts(4, end) = 2; % cell_cell_repulsion_strength (applies to leaders and followers)
    Other_Rescue_Expts(5, end) = 0.05; % cell_cell_repulsion_strength (applies to leaders and followers)
    Other_Rescue_Expts(6:7, 1) = (42:43)'; % scenario
    Other_Rescue_Expts(6, end) = 2; % cell_cell_repulsion_strength (will apply to followers only)
    Other_Rescue_Expts(7, end) = 0.05; % cell_cell_repulsion_strength (will apply to followers only)
    Other_Rescue_Expts(8:12, 1) = 44*ones(5, 1); % scenario
    Other_Rescue_Expts(9, 4) = 10; % average_time_until_next_filopodia_drop (min)
    Other_Rescue_Expts(10, 1) = 45; %scenario
    Other_Rescue_Expts(11, end) = 2; % cell_cell_repulsion_strength (applies to leaders and followers)
    Other_Rescue_Expts(12, end) = 0.05; % cell_cell_repulsion_strength (applies to leaders and followers)
    Other_Rescue_Expts(13:14, 1) = (46:47)'; % scenario
    Other_Rescue_Expts(13, end) = 2; % cell_cell_repulsion_strength (will apply to followers only)
    Other_Rescue_Expts(14, end) = 0.05; % cell_cell_repulsion_strength (will apply to followers only)
    Other_Rescue_Expts(15, 1) = 31; % scenario
    Other_Rescue_Expts(15, 4) = 10; % average_time_until_next_filopodia_drop (min)
    Other_Rescue_Expts(16, 1) = 48; % scenario
    Other_Rescue_Expts(17:19, 3) = 0.33*ones(3,1); % bias_toward_VM_direction
    Other_Rescue_Expts(17, 4) = 10; % average_time_until_next_filopodia_drop (min)
    Other_Rescue_Expts(17:19, end) = zeros(3,1); % cell_cell_repulsion_strength (applies to leaders and followers)
    Other_Rescue_Expts(18, 1) = 29; % scenario
    Other_Rescue_Expts(19, 1) = 31; % scenario
    Other_Rescue_Expts(20, 3) = 0.25; % bias_toward_VM_direction
    Other_Rescue_Expts(20, 4) = 10; % average_time_until_next_filopodia_drop (min)
    
    Other_Rescue_Expts = repmat(Other_Rescue_Expts, 3, 1);
    Other_Rescue_Expts(21:40, 2) = 60*ones(20, 1); % fibronectin_puncta_wavelength
    
    Other_Rescue_Expts(41:60, 1) = 20*ones(20, 1); % scenario
    Other_Rescue_Expts(41:45,1) = 49*ones(5,1); % scenario
    Other_Rescue_Expts(42, 1) = 58; % scenario
    Other_Rescue_Expts(42, 4) = 10; % average_time_until_next_filopodia_drop (min)
    Other_Rescue_Expts(43, 1) = 50; %scenario
    Other_Rescue_Expts(44, end) = 2; % cell_cell_repulsion_strength (applies to leaders and followers)
    Other_Rescue_Expts(45, end) = 0.05; % cell_cell_repulsion_strength (applies to leaders and followers)
    Other_Rescue_Expts(46:47, 1) = (51:52)'; % scenario
    Other_Rescue_Expts(46, end) = 2; % cell_cell_repulsion_strength (will apply to followers only)
    Other_Rescue_Expts(47, end) = 0.05; % cell_cell_repulsion_strength (will apply to followers only)
    Other_Rescue_Expts(48:52, 1) = 53*ones(5, 1); % scenario
    Other_Rescue_Expts(49, 1) = 59; % scenario
    Other_Rescue_Expts(49, 4) = 10; % average_time_until_next_filopodia_drop (min)
    Other_Rescue_Expts(50, 1) = 54; %scenario
    Other_Rescue_Expts(51, end) = 2; % cell_cell_repulsion_strength (applies to leaders and followers)
    Other_Rescue_Expts(52, end) = 0.05; % cell_cell_repulsion_strength (applies to leaders and followers)
    Other_Rescue_Expts(53:54, 1) = (55:56)'; % scenario
    Other_Rescue_Expts(53, end) = 2; % cell_cell_repulsion_strength (will apply to followers only)
    Other_Rescue_Expts(54, end) = 0.05; % cell_cell_repulsion_strength (will apply to followers only)
    Other_Rescue_Expts(55, 1) = 60; % scenario
    Other_Rescue_Expts(55, 4) = 10; % average_time_until_next_filopodia_drop (min)
    Other_Rescue_Expts(56, 1) = 57; % scenario
    Other_Rescue_Expts(57:59, 3) = 0.33*ones(3,1); % bias_toward_VM_direction
    Other_Rescue_Expts(57, 1) = 33; % scenario
    Other_Rescue_Expts(57, 4) = 10; % average_time_until_next_filopodia_drop (min)
    Other_Rescue_Expts(57:59, end) = zeros(3,1); % cell_cell_repulsion_strength (applies to leaders and followers)
    Other_Rescue_Expts(58, 1) = 34; % scenario
    Other_Rescue_Expts(59, 1) = 31; % scenario
    Other_Rescue_Expts(60, 1) = 33; % scenario
    Other_Rescue_Expts(60, 3) = 0.25; % bias_toward_VM_direction
    Other_Rescue_Expts(60, 4) = 10; % average_time_until_next_filopodia_drop (min)
    
    parameters.values = [parameters.values; Other_Rescue_Expts];
    
    Follow_Up_Expts = repmat(parameters.values(1,:), 10, 1);
    Follow_Up_Expts(1, 4) = 5; % average_time_until_next_filopodia_drop (min)
    Follow_Up_Expts(2, 4) = 1; % average_time_until_next_filopodia_drop (min)
    Follow_Up_Expts(3:6, 3) = ones(4,1); % bias_toward_VM_direction
    Follow_Up_Expts(4, 4) = 10; % average_time_until_next_filopodia_drop (min)
    Follow_Up_Expts(5:6, 1) = 17*ones(2,1); % scenario
    Follow_Up_Expts(6, 4) = 10; % average_time_until_next_filopodia_drop (min)
    Follow_Up_Expts(7, 1) = 4; % scenario
    Follow_Up_Expts(7:10, 2) = 60*ones(4,1); % fibronectin_puncta_wavelength
    
    parameters.values = [parameters.values; Follow_Up_Expts];
    
    Final_Table_Expts = repmat(parameters.values(1,:), 35, 1);
    Final_Table_Expts(1,1) = 61; %scenario
    Final_Table_Expts(2,1) = 62; %scenario
    Final_Table_Expts(3:4,1) = 63*ones(2,1); %scenario
    Final_Table_Expts(1:4,3) = 0.33*ones(4,1); % bias_toward_VM_direction
    Final_Table_Expts(3,end) = 0.05; % cell_cell_repulsion_strength (will apply to followers only)
    Final_Table_Expts(4,end) = 2; % cell_cell_repulsion_strength (will apply to followers only)
    Final_Table_Expts(5,1) = 64; %scenario
    Final_Table_Expts(6,1) = 65; %scenario
    Final_Table_Expts(7,1) = 66; %scenario
    Final_Table_Expts(8,1) = 66; %scenario
    Final_Table_Expts(7,end) = 0.05; % cell_cell_repulsion_strength (will apply to followers only)
    Final_Table_Expts(8,end) = 2; % cell_cell_repulsion_strength (will apply to followers only)
    Final_Table_Expts(9,1) = 61; %scenario
    Final_Table_Expts(10,1) = 62; %scenario
    Final_Table_Expts(11,1) = 63; %scenario
    Final_Table_Expts(12,1) = 63; %scenario
    Final_Table_Expts(9:12,4) = 10*ones(4,1); % average_time_until_next_filopodia_drop (min)
    Final_Table_Expts(11,end) = 0.05; % cell_cell_repulsion_strength (will apply to followers only)
    Final_Table_Expts(12,end) = 2; % cell_cell_repulsion_strength (will apply to followers only)
    Final_Table_Expts(13,1) = 29; %scenario
    Final_Table_Expts(14,1) = 29; %scenario
    Final_Table_Expts(13,end) = 0.05; % cell_cell_repulsion_strength (will apply to followers only)
    Final_Table_Expts(14,end) = 2; % cell_cell_repulsion_strength (will apply to followers only)
    Final_Table_Expts(15,1) = 67; %scenario
    Final_Table_Expts(16,1) = 68; %scenario
    Final_Table_Expts(17:18,1) = 69*ones(2,1); %scenario
    Final_Table_Expts(17,end) = 0.05; % cell_cell_repulsion_strength (will apply to followers only)
    Final_Table_Expts(18,end) = 2; % cell_cell_repulsion_strength (will apply to followers only)
    Final_Table_Expts(19,1) = 61; %scenario
    Final_Table_Expts(20,1) = 62; %scenario
    Final_Table_Expts(21:22,1) = 63*ones(2,1); %scenario
    Final_Table_Expts(21,end) = 0.05; % cell_cell_repulsion_strength (will apply to followers only)
    Final_Table_Expts(22,end) = 2; % cell_cell_repulsion_strength (will apply to followers only)
    Final_Table_Expts(23,1) = 70; %scenario
    Final_Table_Expts(24:25, 1) = 31*ones(2,1); %scenario
    Final_Table_Expts(24,end) = 0.05; % cell_cell_repulsion_strength (will apply to both cells)
    Final_Table_Expts(25,end) = 2; % cell_cell_repulsion_strength (will apply to both cells)
    Final_Table_Expts(26,1) = 71; %scenario
    Final_Table_Expts(27,1) = 72; %scenario
    Final_Table_Expts(28:29,1) = 73*ones(2,1); %scenario
    Final_Table_Expts(28,end) = 0.05; % cell_cell_repulsion_strength (will apply to followers only)
    Final_Table_Expts(29,end) = 2; % cell_cell_repulsion_strength (will apply to followers only)
    Final_Table_Expts(30,1) = 44; %scenario
    Final_Table_Expts(30,3) = 0.33; % bias_toward_VM_direction
    Final_Table_Expts(31,1) = 74; %scenario
    Final_Table_Expts(32,1) = 75; %scenario
    Final_Table_Expts(33,1) = 76; %scenario
    Final_Table_Expts(34:35,1) = 77*ones(2,1); %scenario
    Final_Table_Expts(34,end) = 0.05; % cell_cell_repulsion_strength (will apply to followers only)
    Final_Table_Expts(35,end) = 2; % cell_cell_repulsion_strength (will apply to followers only)
    
    parameters.values = [parameters.values; Final_Table_Expts];
    
    save(strcat(save_path, 'experiments_parameter_matrix_with_LHS_nov2021.mat'), 'parameters');
else
    load(strcat(save_path, 'experiments_parameter_matrix_with_LHS_nov2021.mat'));
end

number_of_parameter_regimes = size(parameters.values, 1);
number_of_realizations_per_parameter_regime = 200; 
maximum_realization_time_per_realization = 1.5; % min
% To estimate runtime, add on an extra parameter regime to reflect the fact
% that, for the final scenario, you run for 1440 min instead of 720 min:
est_runtime_hours = (number_of_parameter_regimes+1)*(number_of_realizations_per_parameter_regime*maximum_realization_time_per_realization+1)/60; % Add on an extra minute for producing the video files and generating the summary statistics
est_runtime_days = est_runtime_hours/24;

fprintf(strcat('The estimated maximum amount of time that this function will take to complete is \t', num2str(est_runtime_hours),' hours or \t', num2str(est_runtime_days), ' days.\n'));
input('Do you still wish to continue (press Enter if yes, Ctrl+C if no)?');

Values_for_this_parameter_set = cell(1,number_of_parameters);
Other_Values = cell(1, 3);
Chemotaxis_Parameter_Names = {'leaders_have_global_signal', 'followers_have_global_signal', 'weight_towards_chemotaxis_cues'};
%% Loop over each parameter value, edit the XML config file, and run for multiple realizations:
for j = 1:number_of_parameter_regimes
    if isfile(strcat(save_path,num2str(j, '%04.f'),'.tar.gz'))
        continue; % skip if you have already generated this data set
    end
    for col = 1:number_of_parameters
        Values_for_this_parameter_set{col} = parameters.values(j, col);
    end % for col
    Parameters = {parameters.names{1:number_of_parameters}; Values_for_this_parameter_set{1:number_of_parameters}};
    py_file = ' ./edit_PhysiCell_settings_XML.py';
    cmd = strcat('python3 ', py_file, sprintf(' %s %d', Parameters{:}));  
    % Use a python function to edit the config file
    system(cmd);       
    if j > 135 % if You are in a parameter regime in which you're adjusting the global signal parameter:
        if j < 160
            for col = 1:3
                Other_Values{col} = Chemotaxis_Parameters(j-135, col);
            end
        elseif ((j==160) || (j == 182) || (j == 204))
            Other_Values = {0, 0, 0};
        elseif j == 293
            Other_Values = {1, 1, 0.25};
        elseif j == 294
            Other_Values = {1, 1, 0.5};
        elseif j == 295
            Other_Values = {1, 1, 0.75};
        else
            Other_Values = {1, 0, 0.75};
        end
        Other_Parameters = {Chemotaxis_Parameter_Names{1:3}; Other_Values{1:3}};
        py_file = ' ./edit_PhysiCell_settings_XML_for_global_signal.py';
        cmd = strcat('python3 ', py_file, sprintf(' %s %d', Other_Parameters{:}));
        system(cmd);
    end
    fprintf(['Parameter set number ',num2str(j), ' of ',num2str(number_of_parameter_regimes),'...\n']);
    %% Compute all of the summary statistics you require:                         
    Summary_Statistics = NCC_Model_Puncta_Multiple_Realizations_experiments(j, number_of_realizations_per_parameter_regime);
end % for j
