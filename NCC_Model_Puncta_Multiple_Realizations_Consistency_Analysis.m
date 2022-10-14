function Summary_Statistics = NCC_Model_Puncta_Multiple_Realizations_Consistency_Analysis(distribution_number, number_of_realizations)
% If you have a MAC, make sure that you specify the correct version of the
% PATH variable so that you have the right version of python 3.x:
str = computer;
desktop = true;
if str(1) == 'M'
    desktop = false;
    setenv('PATH', '/usr/local/bin:/usr/local/sbin:/usr/bin:/bin:/usr/sbin:/sbin');
end
%% Pre-allocate data structures:
% Assuming x = [0, 500] um, and 20 um histogram partitions:
number_of_domain_partitions = 25;
edges = linspace(0, 500, number_of_domain_partitions+1);
Density_Radii = (5:5:100)'; % 20 radii to test

Sum_of_All_Unit_Velocities_for_Order_Parameter_Calculation = zeros(number_of_realizations, 3);% new as of 8/3/21
Sum_of_Leader_Unit_Velocities_for_Order_Parameter_Calculation = zeros(number_of_realizations, 3);% new as of 8/3/21
Sum_of_Follower_Unit_Velocities_for_Order_Parameter_Calculation = zeros(number_of_realizations, 3); % new as of 8/3/21
Speed_Second_Moment = zeros(number_of_realizations, 1); % new as of 8/3/21
Leader_x_position_Second_Moment = zeros(number_of_realizations, 1); % new as of 8/3/21
Follower_x_position_Second_Moment = zeros(number_of_realizations, 1); % new as of 8/3/21
max_cells = 250; % Assume that no more than 250 cells are created each realization
Min_Distances = zeros(max_cells,1); 
Number_of_Cells_within_50_um_Ball = zeros(max_cells,1);
cell_index = 0;

Summary_Statistics = struct('Number_of_Realizations', number_of_realizations,...
    'Number_of_Leaders', zeros(number_of_realizations, 1),...
    'Number_of_Followers', zeros(number_of_realizations, 1), ...
    'Number_of_Oriented_FN_Puncta', zeros(number_of_realizations, 1), ... % new as of 8/3/21
    'Order_Parameter_All_Cells', zeros(number_of_realizations, 1),... % new as of 8/3/21
    'Order_Parameter_Leaders', zeros(number_of_realizations, 1),... 
    'Order_Parameter_Followers', zeros(number_of_realizations, 1), ...
    'Max_X_Length_of_Stream', zeros(number_of_realizations, 1),... % new as of 8/3/21
    'Max_Y_Length_of_Stream', zeros(number_of_realizations, 1), ...
    'Average_FN_Puncta_Orientations', zeros(number_of_realizations, 1), ... %new as of 8/3/21
    'Standard_Deviation_FN_Puncta_Orientations', zeros(number_of_realizations, 1), ... %new as of 8/3/21
    'Est_Stream_Length', zeros(number_of_realizations, 1), ... %new as of 8/3/21
    'Est_Stream_Width', zeros(number_of_realizations, 1), ... %new as of 8/3/21
    'Average_Speed', zeros(number_of_realizations, 1),...
    'Leader_Cell_Average_x_position', zeros(number_of_realizations, 1),...
    'Follower_Cell_Average_x_position', zeros(number_of_realizations, 1),...
    'Leader_Cell_Standard_Deviation_x_position', zeros(number_of_realizations, 1),...
    'Follower_Cell_Standard_Deviation_x_position', zeros(number_of_realizations, 1),...
    'Average_Min_Distance_per_Realization', zeros(number_of_realizations, 1),... % new as of 8/3/21
    'Standard_Deviation_Min_Distance_per_Realization', zeros(number_of_realizations, 1),... % new as of 8/3/21
    'Average_Number_of_Cells_within_50_um_Ball_per_Realization', zeros(number_of_realizations, 1),... % new as of 8/3/21
    'Std_Dev_Number_of_Cells_within_50_um_Ball_per_Realization', zeros(number_of_realizations, 1),... % new as of 8/3/21
    'Histogram_Leader_Cell_x_Positions', zeros(number_of_domain_partitions, number_of_realizations),...
    'Histogram_Follower_Cell_x_Positions', zeros(number_of_domain_partitions, number_of_realizations),...
    'Histogram_Average_Leader_Cell_x_Position', zeros(number_of_domain_partitions, 1),...
    'Histogram_Standard_Deviation_Leader_Cell_x_Position', zeros(number_of_domain_partitions, 1),...
    'Histogram_Average_Follower_Cell_x_Position', zeros(number_of_domain_partitions, 1),...
    'Histogram_Standard_Deviation_Follower_Cell_x_Position', zeros(number_of_domain_partitions, 1),...
    'Histogram_Minimum_Distance_Measure', zeros(number_of_domain_partitions, number_of_realizations),...
    'Histogram_Average_Minimum_Distance_Measure', zeros(number_of_domain_partitions, 1),...
    'Histogram_Standard_Deviation_Minimum_Distance_Measure', zeros(number_of_domain_partitions, 1),...
    'Density_Radii', Density_Radii,...
    'Histogram_Cells_within_Radius_Statistic', zeros(number_of_domain_partitions, length(Density_Radii), number_of_realizations),...
    'Histogram_Average_Cells_within_Radius_Statistic', zeros(number_of_domain_partitions, length(Density_Radii)),...
    'Histogram_Standard_Deviation_Cells_within_Radius_Statistic', zeros(number_of_domain_partitions, length(Density_Radii)),...
    'Persistence_Homology_Largest_Radius_Before_All_Connected', zeros(number_of_realizations,1), ...
    'Gap_Statistic', zeros(number_of_realizations, 1), ...
    'Gap_Statistic_Standard_1_Error', zeros(number_of_realizations, 1), ...
    'Histogram_Bin_Edges', edges');

%     'Average_Leader_Angle', zeros(number_of_realizations, 1),...
%     'Standard_Deviation_Leader_Angle', zeros(number_of_realizations, 1),...
%     'Average_Follower_Angle', zeros(number_of_realizations, 1),...
%     'Standard_Deviation_Follower_Angle', zeros(number_of_realizations, 1),...
%% Run PhysiCell for as many realizations as you desire:
for i = 1:number_of_realizations
    report = strcat("Number of Realizations Run = ", num2str(i), '/', num2str(number_of_realizations),'...\n');
    fprintf(report);

    % Run PhysiCell to get all relevant data for the current realization:
    system("./project-NCC-puncta-experiment");
    fprintf("\n\n\n");

    %% Record all data in the Summary_Statistics structure.
    listing = dir('./output/*.xml');
    MCDS = read_MultiCellDS_xml(listing(end).name, './output');
    % Collect ID labels:
    Leaders = (MCDS.discrete_cells.metadata.type == 0);
    Followers = (MCDS.discrete_cells.metadata.type == 1);
    leaders_and_followers = logical(Leaders + Followers);
    Puncta = (MCDS.discrete_cells.metadata.type == 2);
    viable_puncta = logical(Puncta.*MCDS.discrete_cells.custom.time_since_last_filopodia_drop);
    
    % Count total number of Cells:
    Summary_Statistics.Number_of_Leaders(i) = sum(Leaders);
    Summary_Statistics.Number_of_Followers(i) = sum(Followers);
    Summary_Statistics.Number_of_Oriented_FN_Puncta(i) = sum(viable_puncta);
    
    % Collect sum of unit velocities for the global order parameter
    % average and standard deviation:
    Sum_of_All_Unit_Velocities_for_Order_Parameter_Calculation(i,:) = sum( MCDS.discrete_cells.custom.total_velocity(leaders_and_followers,:)./vecnorm(MCDS.discrete_cells.custom.total_velocity(leaders_and_followers,:),2,2) );
    Sum_of_Leader_Unit_Velocities_for_Order_Parameter_Calculation(i,:) = sum( MCDS.discrete_cells.custom.total_velocity(Leaders,:)./vecnorm(MCDS.discrete_cells.custom.total_velocity(Leaders,:),2,2) );
    Sum_of_Follower_Unit_Velocities_for_Order_Parameter_Calculation(i,:) = sum( MCDS.discrete_cells.custom.total_velocity(Followers,:)./vecnorm(MCDS.discrete_cells.custom.total_velocity(Followers,:),2,2) );
    
    % Also collect order parameter average for each realization:
    Summary_Statistics.Order_Parameter_All_Cells(i) = vecnorm( mean( MCDS.discrete_cells.custom.total_velocity(leaders_and_followers,:)./vecnorm(MCDS.discrete_cells.custom.total_velocity(leaders_and_followers,:),2,2) ) );
    Summary_Statistics.Order_Parameter_Leaders(i) = vecnorm( mean( MCDS.discrete_cells.custom.total_velocity(Leaders,:)./vecnorm(MCDS.discrete_cells.custom.total_velocity(Leaders,:),2,2) ) );
    if Summary_Statistics.Number_of_Followers(i) > 0
        Summary_Statistics.Order_Parameter_Followers(i) = vecnorm( mean( MCDS.discrete_cells.custom.total_velocity(Followers,:)./vecnorm(MCDS.discrete_cells.custom.total_velocity(Followers,:),2,2) ) );
    end
    % Collect maximum x-length and y-width of the stream (NOTE: NOT the
    % Stream length or width if the orientation of the stream is NOT 0
    % rad!)
    Summary_Statistics.Max_X_Length_of_Stream(i) = max(MCDS.discrete_cells.state.position(leaders_and_followers, 1));
    Summary_Statistics.Max_Y_Length_of_Stream(i) = max(MCDS.discrete_cells.state.position(leaders_and_followers, 2)) - min(MCDS.discrete_cells.state.position(leaders_and_followers, 2));
    
    % Collect SUM of speeds for each realization:
    Summary_Statistics.Average_Speed(i) = sum( vecnorm( MCDS.discrete_cells.custom.total_velocity(leaders_and_followers,:), 2, 2) );
    Speed_Second_Moment(i) = sum( vecnorm( MCDS.discrete_cells.custom.total_velocity(leaders_and_followers,:), 2, 2).^2 );
    
    % Collect SUM of angles for each realization so that you can compute a global metric also:
    if Summary_Statistics.Number_of_Oriented_FN_Puncta(i)>0
        Summary_Statistics.Average_FN_Puncta_Orientations(i) = sum(MCDS.discrete_cells.custom.time_until_next_filopodia_drop(viable_puncta));
        Summary_Statistics.Standard_Deviation_FN_Puncta_Orientations(i) = std(MCDS.discrete_cells.custom.time_until_next_filopodia_drop(viable_puncta));
        est_angle = mean(MCDS.discrete_cells.custom.time_until_next_filopodia_drop(viable_puncta));
    else
        est_angle = 0;
    end
    rot = [cos(est_angle), -sin(est_angle); sin(est_angle), cos(est_angle)];
    uv_plane = rot*(MCDS.discrete_cells.state.position(leaders_and_followers, 1:2)');
    Summary_Statistics.Est_Stream_Length(i) = max(uv_plane(1,:));
    Summary_Statistics.Est_Stream_Width(i) = max(uv_plane(2,:))-min(uv_plane(2,:)); 
    
    % Record SUM of positions before you compute the average (but keep the
    % standard deviations for each realization)
    Summary_Statistics.Leader_Cell_Average_x_position(i) = sum(MCDS.discrete_cells.state.position(Leaders, 1));
    Leader_x_position_Second_Moment(i) = sum(MCDS.discrete_cells.state.position(Leaders, 1).^2);
    Summary_Statistics.Follower_Cell_Average_x_position(i) = sum(MCDS.discrete_cells.state.position(Followers, 1));
    Follower_x_position_Second_Moment(i) = sum(MCDS.discrete_cells.state.position(Followers, 1).^2);
    Summary_Statistics.Leader_Cell_Standard_Deviation_x_position(i) = std(MCDS.discrete_cells.state.position(Leaders, 1));
    if Summary_Statistics.Number_of_Followers(i) > 0
        Summary_Statistics.Follower_Cell_Standard_Deviation_x_position(i) = std(MCDS.discrete_cells.state.position(Followers, 1));
    end
    Summary_Statistics.Histogram_Leader_Cell_x_Positions(:, i) = histcounts(MCDS.discrete_cells.state.position(Leaders, 1), edges)';
    Summary_Statistics.Histogram_Follower_Cell_x_Positions(:, i) = histcounts(MCDS.discrete_cells.state.position(Followers, 1), edges)';
    
    
    % Reduce the position vector to only include leaders and followers:
    MCDS.discrete_cells.state.position = MCDS.discrete_cells.state.position(leaders_and_followers, :);
    
    % Get a count of all of the cells in each histogram bin:
    cell_numbers = histcounts(MCDS.discrete_cells.state.position(:, 1), edges)';
    if sum(cell_numbers) > max_cells
        max_cells = sum(cell_numbers);
        Min_Distances = zeros(max_cells,1); 
        Number_of_Cells_within_50_um_Ball = zeros(max_cells,1);
    end
    
    for ii = 1:length(cell_numbers)
        % Get all cell indices within the current histogram bin
        cell_indices_within_bin = find( (MCDS.discrete_cells.state.position(:,1) >= edges(ii)) == (MCDS.discrete_cells.state.position(:,1) < edges(ii+1)) );
        % If there are no cells, statistics are set equal to 0:
        if isempty(cell_indices_within_bin) || length(cell_indices_within_bin) == 1
           Summary_Statistics.Histogram_Minimum_Distance_Measure(ii, i) = 0; 
           Summary_Statistics.Histogram_Cells_within_Radius_Statistic(ii, :, i) = 0;  
           continue; 
        end
        % Create a vector for the minimum distance between cells (value should be
        % larger if the stream is more disperse)
        min_distance_for_each_cell = zeros(length(cell_indices_within_bin), 1);
        number_of_cells_for_each_persistence_radius = zeros(length(cell_indices_within_bin), length(Density_Radii));
        for j = 1:length(cell_indices_within_bin)
            % For each cell, calculate the shortest distance with any other
            % cell. Make sure to exclude yourself.
            distances_for_all_cells_but_myself = sqrt( sum( ([MCDS.discrete_cells.state.position(1:cell_indices_within_bin(j)-1,:); MCDS.discrete_cells.state.position(cell_indices_within_bin(j)+1:end,:)]- MCDS.discrete_cells.state.position(cell_indices_within_bin(j),:) ).^2, 2));
            min_distance_for_each_cell(j) = min(distances_for_all_cells_but_myself);
            % For each cell and persistence homology radius "eps", count 
            % the number of cells that appear in a ball of radius eps 
            % around cell j:
            for r = 1:length(Density_Radii)
                number_of_cells_for_each_persistence_radius(j, r) = sum((distances_for_all_cells_but_myself<=Density_Radii(r)));
            end % for r
        end % for j
        Min_Distances(cell_index+1:cell_index+length(min_distance_for_each_cell)) = min_distance_for_each_cell(:)+1e-31;
        Number_of_Cells_within_50_um_Ball(cell_index+1:cell_index+length(min_distance_for_each_cell)) = number_of_cells_for_each_persistence_radius(:, 10) + 1e-31;
        cell_index = cell_index + length(min_distance_for_each_cell);
        Summary_Statistics.Histogram_Minimum_Distance_Measure(ii,i) = mean(min_distance_for_each_cell);
        Summary_Statistics.Histogram_Cells_within_Radius_Statistic(ii, :, i) = mean(number_of_cells_for_each_persistence_radius, 1);
%         if isnan(Summary_Statistics.Histogram_Minimum_Distance_Measure(ii,i))
%             Summary_Statistics.Histogram_Minimum_Distance_Measure(ii,i) = 0;
%         end
%         if isnan(Summary_Statistics.Histogram_Cells_within_Radius_Statistic(ii, :, i))
%            Summary_Statistics.Histogram_Cells_within_Radius_Statistic(ii, :, i) = 0; 
%         end
    end
    Summary_Statistics.Average_Min_Distance_per_Realization(i) = mean(Min_Distances(Min_Distances>0)-1e-31);
    Summary_Statistics.Standard_Deviation_Min_Distance_per_Realization(i) = std(Min_Distances(Min_Distances>0)-1e-31);
    Summary_Statistics.Average_Number_of_Cells_within_50_um_Ball_per_Realization(i) = mean(Number_of_Cells_within_50_um_Ball(Number_of_Cells_within_50_um_Ball>0)-1e-31);
    Summary_Statistics.Std_Dev_Number_of_Cells_within_50_um_Ball_per_Realization(i) = std(Number_of_Cells_within_50_um_Ball(Number_of_Cells_within_50_um_Ball>0)-1e-31);
    Min_Distances(Min_Distances>0) = 0;
    Number_of_Cells_within_50_um_Ball(Number_of_Cells_within_50_um_Ball>0) = 0;
    cell_index = 0;
    % NOTICE: Remember that you have reduced the position vector of the
    %           MCDS structure to ONLY include leaders and followers!
%     % Write a csv file with the pairwise distances between cells:
%     dist_matrix = squareform(pdist(MCDS.discrete_cells.state.position));
    try
        writematrix(MCDS.discrete_cells.state.position, 'final_cell_positions_puncta_model.csv');
    catch
        csvwrite('final_cell_positions_puncta_model.csv', MCDS.discrete_cells.state.position);
    end
%     writematrix(dist_matrix, 'pairwise_cell_distances.csv');
    % Compute ripser diagrams and load them into matlab:
    system('python3 compute_persistence_homology_diagrams.py');
    dgms = csvread('persistence_homology_diagrams_puncta_model.csv');
    % WARNING: If you're computing using Vietoris-Rips complex, the output
    % is the DIAMETERS, NOT THE RADII, of the balls used to compute the
    % filtration. Divide by 2 to correct for this.
    
    % Store the next to last RADIUS of the connected components:
    Summary_Statistics.Persistence_Homology_Largest_Radius_Before_All_Connected(i) = dgms(end-1,2)./2;
    % Compute the gap statistic, using the 1-standard error method:
    system('python3 compute_gap_statistic.py');
    Summary_Statistics.Gap_Statistic(i) = csvread('gap_statistic_puncta_model.txt');
    Summary_Statistics.Gap_Statistic_Standard_1_Error(i) = csvread('gap_std_1_error_statistic_puncta_model.txt');
    fprintf("\n\n\n");
end
%% Record all average/standard deviation data over multiple realizations:
Summary_Statistics.Average_FN_Puncta_Orientations = Summary_Statistics.Average_FN_Puncta_Orientations./Summary_Statistics.Number_of_Oriented_FN_Puncta;
Summary_Statistics.Average_FN_Puncta_Orientations(isnan(Summary_Statistics.Average_FN_Puncta_Orientations)) = 0;
Summary_Statistics.Average_FN_Puncta_Orientations(isinf(Summary_Statistics.Average_FN_Puncta_Orientations)) = 0;

Summary_Statistics.Average_Speed = Summary_Statistics.Average_Speed./(Summary_Statistics.Number_of_Leaders + Summary_Statistics.Number_of_Followers);

% Averages for histogram data -- averages are over all of the bins, not the
% realizations
Summary_Statistics.Histogram_Average_Leader_Cell_x_Position = mean(Summary_Statistics.Histogram_Leader_Cell_x_Positions, 2);
Summary_Statistics.Histogram_Standard_Deviation_Leader_Cell_x_Position = std(Summary_Statistics.Histogram_Leader_Cell_x_Positions, 0, 2);
Summary_Statistics.Histogram_Average_Follower_Cell_x_Position = mean(Summary_Statistics.Histogram_Follower_Cell_x_Positions, 2);
Summary_Statistics.Histogram_Standard_Deviation_Follower_Cell_x_Position = std(Summary_Statistics.Histogram_Follower_Cell_x_Positions, 0, 2);

% The following two lines represent an average/standard deviation of 
% averages:
Summary_Statistics.Histogram_Average_Minimum_Distance_Measure = mean(Summary_Statistics.Histogram_Minimum_Distance_Measure, 2);
Summary_Statistics.Histogram_Standard_Deviation_Minimum_Distance_Measure = std(Summary_Statistics.Histogram_Minimum_Distance_Measure, 0, 2);

% The following two lines represent the average/standard deviation of
% averages:
% Summary_Statistics.Histogram_Average_Cells_within_Radius_Statistic = mean(Summary_Statistics.Histogram_Cells_within_Radius_Statistic, 3);
% Summary_Statistics.Histogram_Standard_Deviation_Cells_within_Radius_Statistic = std(Summary_Statistics.Histogram_Cells_within_Radius_Statistic, 0, 3);

% Since we know how many leaders and followers are in each x-bin
% for every realization, we can get rid of the inner averages to get a 
% proper statistic for the average number of cells within a 50 um ball, put into a histogram:
cells_in_every_bin_per_realization = permute(repmat(Summary_Statistics.Histogram_Leader_Cell_x_Positions + Summary_Statistics.Histogram_Follower_Cell_x_Positions, 1, 1, length(Density_Radii)), [1 3 2]);
% Total number of cell counts for all realizations, divided by number of
% cells in each realization:
Summary_Statistics.Histogram_Average_Cells_within_Radius_Statistic = sum(Summary_Statistics.Histogram_Cells_within_Radius_Statistic.*cells_in_every_bin_per_realization, 3)./sum(cells_in_every_bin_per_realization, 3);
% Set all NaN values to 0:
Summary_Statistics.Histogram_Average_Cells_within_Radius_Statistic(isnan(Summary_Statistics.Histogram_Average_Cells_within_Radius_Statistic)) = 0;
% Standard deviation of all the cell counts across all realizations:
Summary_Statistics.Histogram_Standard_Deviation_Cells_within_Radius_Statistic = sqrt(1./(sum(cells_in_every_bin_per_realization, 3)-1)).*sqrt(sum((Summary_Statistics.Histogram_Cells_within_Radius_Statistic.*cells_in_every_bin_per_realization - Summary_Statistics.Histogram_Average_Cells_within_Radius_Statistic).^2, 3));
% Get rid of NaN, Inf Values:
Summary_Statistics.Histogram_Standard_Deviation_Cells_within_Radius_Statistic(isnan(Summary_Statistics.Histogram_Standard_Deviation_Cells_within_Radius_Statistic)) = 0;
Summary_Statistics.Histogram_Standard_Deviation_Cells_within_Radius_Statistic(isinf(Summary_Statistics.Histogram_Standard_Deviation_Cells_within_Radius_Statistic)) = 0;

%% Save the statistics in the output folder, and create a csv file for the most relevant summary statistics:
% NOTE: Currently we only store the data for the last realization, not the
% ones in between! 
    % You already took into account the sum in the numerator for the
    % following fraction:
    Summary_Statistics.Leader_Cell_Average_x_position = Summary_Statistics.Leader_Cell_Average_x_position./Summary_Statistics.Number_of_Leaders;
    if sum(Summary_Statistics.Number_of_Followers)>0
        % You already took into account the sum in the numerator for the
        % following fraction:
        Summary_Statistics.Follower_Cell_Average_x_position = Summary_Statistics.Follower_Cell_Average_x_position./Summary_Statistics.Number_of_Followers;
    end
video_title = strcat(num2str(number_of_realizations),'_',num2str(distribution_number));
save(['./output/',video_title,'_summary_statistics.mat'], 'Summary_Statistics');
    %% Transfer Everything to a New Directory and compress the resulting folder:
    system(strcat('mkdir ./NCC_consistency_analysis/', video_title)); 
    system(strcat('mv ./output/* ./NCC_consistency_analysis/',video_title) );
    system('touch ./output/empty.txt');
    system(['cd ./NCC_consistency_analysis/',video_title,' && tar -zcf ',video_title,'.tar.gz * && mv ',video_title,'.tar.gz ../ && cd -']);
    system(['rm -rf ./NCC_consistency_analysis/',video_title]);
end