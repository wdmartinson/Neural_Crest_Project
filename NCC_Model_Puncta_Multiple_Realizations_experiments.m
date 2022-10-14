function Summary_Statistics = NCC_Model_Puncta_Multiple_Realizations_experiments(parameter_set_number, number_of_realizations)
% If you have a MAC, make sure that you specify the correct version of the
% PATH variable so that you have the right version of python 3.x:
desktop = true;
str = computer;
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
max_cells = 250; % Assume that no more than 250 cells are cpreated each realization
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
    'Porosity_Measure', zeros(number_of_realizations, 1), ...
    'Percent_Target_Area_Covered', zeros(number_of_realizations, 1), ...
    'Percent_Cell_Area_Outside_Target_Sites', zeros(number_of_realizations, 1), ...
    'Average_Min_Distance_from_Follower_to_Leader_per_Realization', zeros(number_of_realizations, 1), ... % CAREFUL! THIS STATISTIC IS MISNAMED - IT ACTUALLY MEASURES THE LEADER -> FOLLOWER DISTANCE
    'Average_Min_Distance_from_Leader_to_Follower_per_Realization', zeros(number_of_realizations, 1), ... % CAREFUL! THIS STATISTIC IS MISNAMED - IT ACTUALLY MEASURES THE FOLLOWER -> LEADER DISTANCE
    'Average_Min_Distance_from_Leader_to_Leader_per_Realization', zeros(number_of_realizations, 1), ...
    'Number_of_Single_Leader_Streams', zeros(number_of_realizations, 1), ...
    'Number_of_Single_Leader_Streams_Filtered_by_Density_Stat', zeros(number_of_realizations, 1), ...
    'Number_of_Outlier_Cells', zeros(number_of_realizations, 1),...
    'Max_Leader_X_Position', zeros(number_of_realizations, 1),...
    'Max_Follower_X_Position', zeros(number_of_realizations, 1),...
    'Global_Average_Leader_Speed', zeros(1,1),...
    'Average_Leader_Speed_per_Realization', zeros(number_of_realizations, 1),...
    'Global_Average_Follower_Speed', zeros(1,1),...
    'Average_Follower_Speed_per_Realization', zeros(number_of_realizations, 1),...
    'Min_Leader_X_Position', zeros(number_of_realizations, 1),...
    'Ave_Horizontal_Leader_Speed_per_Realization', zeros(number_of_realizations, 1),...
    'Std_Horizontal_Leader_Speed_per_Realization', zeros(number_of_realizations, 1),...
    'Ave_Horizontal_Follower_Speed_per_Realization', zeros(number_of_realizations, 1),...
    'Std_Horizontal_Follower_Speed_per_Realization', zeros(number_of_realizations, 1),...
    'Ave_Vertical_Leader_Speed_per_Realization', zeros(number_of_realizations, 1),...
    'Std_Vertical_Leader_Speed_per_Realization', zeros(number_of_realizations, 1),...
    'Ave_Vertical_Follower_Speed_per_Realization', zeros(number_of_realizations, 1),...
    'Std_Vertical_Follower_Speed_per_Realization', zeros(number_of_realizations, 1),...
    'Histogram_Bin_Edges', edges');

Xmin = readstruct('./config/PhysiCell_settings.xml').domain.x_min;
Xmax = readstruct('./config/PhysiCell_settings.xml').domain.x_max;
Ymin = readstruct('./config/PhysiCell_settings.xml').domain.y_min;
Ymax = readstruct('./config/PhysiCell_settings.xml').domain.y_max;
l_entr = readstruct('./config/PhysiCell_settings.xml').user_parameters.y_length_of_cell_entrance_strip.Text;
R_filo = readstruct('./config/PhysiCell_settings.xml').user_parameters.filopodia_sensing_radius.Text;

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
    temp = MCDS.discrete_cells.custom.total_velocity(leaders_and_followers,:)./vecnorm(MCDS.discrete_cells.custom.total_velocity(leaders_and_followers,:),2,2);
    temp(isnan(temp)) = 0; % If you had zero velocity, reset those values
    Summary_Statistics.Order_Parameter_All_Cells(i) = vecnorm( mean( temp ) );
    if Summary_Statistics.Number_of_Leaders(i) > 0
        temp = MCDS.discrete_cells.custom.total_velocity(Leaders,:)./vecnorm(MCDS.discrete_cells.custom.total_velocity(Leaders,:),2,2);
        temp(isnan(temp)) = 0;
        Summary_Statistics.Order_Parameter_Leaders(i) = vecnorm( mean( temp ) );
    end
    if Summary_Statistics.Number_of_Followers(i) > 0
        temp = MCDS.discrete_cells.custom.total_velocity(Followers,:)./vecnorm(MCDS.discrete_cells.custom.total_velocity(Followers,:),2,2);
        temp(isnan(temp)) = 0;
        Summary_Statistics.Order_Parameter_Followers(i) = vecnorm( mean( temp ) );
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
    
    % Also collect the average speed of leaders and followers, and their sum for the global statistic:
    Summary_Statistics.Average_Leader_Speed_per_Realization(i) = mean(vecnorm( MCDS.discrete_cells.custom.total_velocity(Leaders,:), 2, 2) );
    Summary_Statistics.Average_Follower_Speed_per_Realization(i) = mean(vecnorm( MCDS.discrete_cells.custom.total_velocity(Followers,:), 2, 2) );
    Summary_Statistics.Global_Average_Leader_Speed = Summary_Statistics.Global_Average_Leader_Speed + sum(vecnorm( MCDS.discrete_cells.custom.total_velocity(Leaders,:), 2, 2) );
    Summary_Statistics.Global_Average_Follower_Speed = Summary_Statistics.Global_Average_Follower_Speed + sum(vecnorm( MCDS.discrete_cells.custom.total_velocity(Followers,:), 2, 2) );
    
    % Also collect the average horizontal and vertical components of the
    % follower and leader velocities:
    Leader_Angles = atan2(MCDS.discrete_cells.custom.total_velocity(Leaders,2), MCDS.discrete_cells.custom.total_velocity(Leaders,1));
    Leader_Speeds = vecnorm( MCDS.discrete_cells.custom.total_velocity(Leaders,:), 2, 2);
    Follower_Angles = atan2(MCDS.discrete_cells.custom.total_velocity(Followers,2), MCDS.discrete_cells.custom.total_velocity(Followers,1));
    Follower_Speeds = vecnorm( MCDS.discrete_cells.custom.total_velocity(Followers,:), 2, 2);
    Summary_Statistics.Ave_Horizontal_Leader_Speed_per_Realization(i) = mean(Leader_Speeds.*cos(Leader_Angles));
    Summary_Statistics.Std_Horizontal_Leader_Speed_per_Realization(i) = std(Leader_Speeds.*cos(Leader_Angles));
    Summary_Statistics.Ave_Vertical_Leader_Speed_per_Realization(i) = mean( Leader_Speeds.*sin(Leader_Angles) );
    Summary_Statistics.Std_Vertical_Leader_Speed_per_Realization(i) = std( Leader_Speeds.*sin(Leader_Angles) );
    Summary_Statistics.Ave_Horizontal_Follower_Speed_per_Realization(i) = mean( Follower_Speeds.*cos(Follower_Angles) );
    Summary_Statistics.Std_Horizontal_Follower_Speed_per_Realization(i) = std( Follower_Speeds.*cos(Follower_Angles) );
    Summary_Statistics.Ave_Vertical_Follower_Speed_per_Realization(i) = mean( Follower_Speeds.*sin(Follower_Angles) );
    Summary_Statistics.Std_Vertical_Follower_Speed_per_Realization(i) = std( Follower_Speeds.*sin(Follower_Angles) );
    
    % Record SUM of positions before you compute the average (but keep the
    % standard deviations for each realization)
    Summary_Statistics.Max_Leader_X_Position(i) = max(MCDS.discrete_cells.state.position(Leaders, 1));
    Summary_Statistics.Min_Leader_X_Position(i) = min(MCDS.discrete_cells.state.position(Leaders, 1));
    Summary_Statistics.Leader_Cell_Average_x_position(i) = sum(MCDS.discrete_cells.state.position(Leaders, 1));
    Leader_x_position_Second_Moment(i) = sum(MCDS.discrete_cells.state.position(Leaders, 1).^2);
    Summary_Statistics.Follower_Cell_Average_x_position(i) = sum(MCDS.discrete_cells.state.position(Followers, 1));
    Follower_x_position_Second_Moment(i) = sum(MCDS.discrete_cells.state.position(Followers, 1).^2);
    Summary_Statistics.Leader_Cell_Standard_Deviation_x_position(i) = std(MCDS.discrete_cells.state.position(Leaders, 1));
    if Summary_Statistics.Number_of_Followers(i) > 0
        Summary_Statistics.Follower_Cell_Standard_Deviation_x_position(i) = std(MCDS.discrete_cells.state.position(Followers, 1));
        Summary_Statistics.Max_Follower_X_Position(i) = max(MCDS.discrete_cells.state.position(Followers, 1));
    end
    Summary_Statistics.Histogram_Leader_Cell_x_Positions(:, i) = histcounts(MCDS.discrete_cells.state.position(Leaders, 1), edges)';
    Summary_Statistics.Histogram_Follower_Cell_x_Positions(:, i) = histcounts(MCDS.discrete_cells.state.position(Followers, 1), edges)';
    
    % Get the porosity measures and use imread to compute the percentage
    % statistics:
    Summary_Statistics.Porosity_Measure(i) = Summary_Statistics.Max_Y_Length_of_Stream(i)./(sum(double(Leaders)) + sum(double(Followers)));
    
    fig3 = figure('units', 'inches', 'position', [0,0,7,5]);
    xlabel('$x$, $\mu$m', 'interpreter', 'latex', 'fontsize', 16);
    ylabel('$y$, $\mu$m', 'interpreter', 'latex', 'fontsize', 16);
    axis square; box on;
    set(gca, 'nextplot', 'replacechildren', 'visible','on');
    ax3 = gca;
    title('','interpreter', 'latex', 'fontsize', 20);
    set(get(ax3, 'title'), 'string', ['$t$ = ',num2str(MCDS.metadata.current_time), ' ', MCDS.metadata.time_units]);
    ax3.Title.Interpreter = 'latex';
    ax3.FontSize = 16;
    xlim([Xmin, Xmax]);
    ylim([Ymin, Ymax]);
    rectangle('Position', [0 -l_entr/2 500 l_entr], 'FaceColor', 'r', 'EdgeColor', 'r');
    saveas(fig3, './output/Rect.png');
    rgb = imread('./output/Rect.png');
    A_box_rect = (rgb(:,:,1) == 255 & rgb(:,:,2) == 0 & rgb(:,:,3) == 0);
    Interest = (abs(MCDS.discrete_cells.state.position(:,2,:)) <= l_entr/2) & leaders_and_followers';
    hold on; interest_radii = (MCDS.discrete_cells.phenotype.geometrical_properties.volumes.nuclear(Interest)./4./pi.*3).^(1/3);
    h1 = viscircles([MCDS.discrete_cells.state.position(Interest,1) MCDS.discrete_cells.state.position(Interest,2)], interest_radii, 'color', 'b', 'linewidth', 0.75);
    if ~isempty(MCDS.discrete_cells.state.position(Leaders,1))
    x_leaders = h1.Children(1).XData;
    y_leaders = h1.Children(1).YData;
    idx = [0, find(isnan(x_leaders))];
    for ii = 1:length(idx)-1
    fill(x_leaders(idx(ii)+1:idx(ii+1)-1), y_leaders(idx(ii)+1:idx(ii+1)-1), 'b', 'facealpha', 1);
    end
    end
    saveas(fig3, './output/Box_Cell.png');
    rgb = imread('./output/Box_Cell.png');
    A_box_cell = (rgb(:,:,1) == 0 & rgb(:,:,2) == 0 & rgb(:,:,3) == 255);
    cell_radii = (MCDS.discrete_cells.phenotype.geometrical_properties.volumes.nuclear(leaders_and_followers)./4./pi.*3).^(1/3);
    h1 = viscircles([MCDS.discrete_cells.state.position(leaders_and_followers,1) MCDS.discrete_cells.state.position(leaders_and_followers,2)], cell_radii, 'color', 'b', 'linewidth', 0.75);
    if ~isempty(MCDS.discrete_cells.state.position(leaders_and_followers,1))
        x_leaders = h1.Children(1).XData;
        y_leaders = h1.Children(1).YData;
        idx = [0, find(isnan(x_leaders))];
        for ii = 1:length(idx)-1
            fill(x_leaders(idx(ii)+1:idx(ii+1)-1), y_leaders(idx(ii)+1:idx(ii+1)-1), 'b', 'facealpha', 1);
        end
    end
    saveas(fig3, './output/Tot_Cell.png');
    rgb = imread('./output/Tot_Cell.png');
    A_tot_cell = (rgb(:,:,1) == 0 & rgb(:,:,2) == 0 & rgb(:,:,3) == 255);
    close all;
    
    Summary_Statistics.Percent_Target_Area_Covered(i) = sum(sum(double(A_box_cell)))./sum(sum(double(A_box_rect)));
    Summary_Statistics.Percent_Cell_Area_Outside_Target_Sites(i) = 1-sum(sum(double(A_box_cell)))./sum(sum(double(A_tot_cell)));
    
    % Record the number of leaders and followers that are outside the
    % target corridor:
    if sum(leaders_and_followers) > 0
        temp_mat = MCDS.discrete_cells.state.position(leaders_and_followers, 2);
        Summary_Statistics.Number_of_Outlier_Cells(i) = sum(abs(temp_mat)>l_entr/2);
    end
    % Record the nearest distance between a leader cell and a follower
    % cell, or between a 
    if ~isempty(MCDS.discrete_cells.state.position(Leaders,1)) && ~isempty(MCDS.discrete_cells.state.position(Followers,1))
        Summary_Statistics.Average_Min_Distance_from_Follower_to_Leader_per_Realization(i) = mean( min( pdist2(MCDS.discrete_cells.state.position(Leaders,:), MCDS.discrete_cells.state.position(Followers,:), 'euclidean') , [], 2) );
        Summary_Statistics.Average_Min_Distance_from_Leader_to_Follower_per_Realization(i) = mean( min( pdist2(MCDS.discrete_cells.state.position(Leaders,:), MCDS.discrete_cells.state.position(Followers,:), 'euclidean') , [], 1) );
    else
        Summary_Statistics.Average_Min_Distance_from_Follower_to_Leader_per_Realization(i) = NaN;
        Summary_Statistics.Average_Min_Distance_from_Leader_to_Follower_per_Realization(i) = NaN;
    end
    if ~isempty(MCDS.discrete_cells.state.position(Leaders,1))
        leader_dist_mat = pdist2(MCDS.discrete_cells.state.position(Leaders,:), MCDS.discrete_cells.state.position(Leaders,:), 'euclidean');
        Summary_Statistics.Average_Min_Distance_from_Leader_to_Leader_per_Realization(i) = mean( min( leader_dist_mat + diag(Inf.*ones(size(leader_dist_mat,1), 1)) , [], 2) );
    else
        Summary_Statistics.Average_Min_Distance_from_Leader_to_Leader_per_Realization(i) = NaN;
    end
    
    
    if Summary_Statistics.Number_of_Followers(i) > 0 && Summary_Statistics.Number_of_Leaders(i) > 0
        D = pdist2(MCDS.discrete_cells.state.position(leaders_and_followers, :), MCDS.discrete_cells.state.position(leaders_and_followers, :), 'euclidean');
        D = D + diag(Inf.*ones(size(D, 1), 1));
        ave_no_Cells_in_50_um_ball = mean([sum((D < 50))'; sum((D<50),2)]);
        % Compute pairwise distance matrix between all cells/puncta
        D = pdist2(MCDS.discrete_cells.state.position, MCDS.discrete_cells.state.position, 'euclidean');
        D = D + diag(Inf.*ones(size(D,1),1));
        % Determine which followers don't have many neighbors:
        Cells_in_50_um_ball_around_followers = sum((D<50).*(Followers & Followers'),2);
        Sparse_Followers = (Cells_in_50_um_ball_around_followers~=0 & Cells_in_50_um_ball_around_followers<ave_no_Cells_in_50_um_ball);
        % Compute the number of outliers (not filtered):
        temp = D.*(Followers'.*viable_puncta);
        idx = (temp~=0 & temp<R_filo); 
        unique_follower_count = 0;
        unique_IDs = [];
        for rr = 1:size(idx,1)
            if ~Followers(rr)
                continue;
            end
            IDs = round(MCDS.discrete_cells.custom.fibronectin_ID(idx(rr,:)));
            if length(unique(IDs(IDs >= 0)))== 1
                unique_follower_count = unique_follower_count + 1;
                unique_IDs = [unique_IDs, unique(IDs(IDs >= 0))];
            end
        end
        [~, idx_IDs, ~] = unique(unique_IDs);
        unique_IDs(idx_IDs) = [];
        Summary_Statistics.Number_of_Single_Leader_Streams(i) = length(unique(unique_IDs));
        % Compute the number of outliers (filtering out dense followers):
        temp = D.*(Sparse_Followers.*viable_puncta);
        idx = (temp~=0 & temp<R_filo);
        filtered_unique_follower_count = 0;
        filtered_unique_IDs = [];
        for rr = 1:size(idx,1)
            IDs = round(MCDS.discrete_cells.custom.fibronectin_ID(idx(rr,:)));
            if length(unique(IDs(IDs >= 0))) == 1
            Leader_pos = MCDS.discrete_cells.state.position(Leaders, :);
            [~, leader_id] = min(pdist2(MCDS.discrete_cells.state.position(rr, :), Leader_pos, 'euclidean'));
            T = pdist2(Leader_pos(leader_id, :), Leader_pos, 'euclidean');
            T(leader_id) = Inf;
            if min(T) <= (R_filo)
                continue;
            end
            filtered_unique_follower_count = filtered_unique_follower_count + 1;
            filtered_unique_IDs = [filtered_unique_IDs, unique(IDs(IDs >= 0))];
            end
        end
        [~, idx_IDs, ~] = unique(filtered_unique_IDs);
        filtered_unique_IDs(idx_IDs) = [];
        Summary_Statistics.Number_of_Single_Leader_Streams_Filtered_by_Density_Stat(i) = length(unique(filtered_unique_IDs));
    end
    
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
    try
        writematrix(MCDS.discrete_cells.state.position, 'final_cell_positions_puncta_model.csv');
    catch
        csvwrite('final_cell_positions_puncta_model.csv', MCDS.discrete_cells.state.position);
    end
    % Compute ripser diagrams and load them into matlab:
    system('python3 compute_persistence_homology_diagrams.py');
    dgms = csvread('persistence_homology_diagrams_puncta_model.csv');
    % WARNING: If you're computing using Vietoris-Rips complex, the output
    % is the DIAMETERS, NOT THE RADII, of the balls used to compute the
    % filtration. Divide by 2 to correct for this.
    % Store the next to last RADIUS of the connected components:
    if size(dgms, 1) > 1
        Summary_Statistics.Persistence_Homology_Largest_Radius_Before_All_Connected(i) = dgms(end-1,2)./2;
    else
        Summary_Statistics.Persistence_Homology_Largest_Radius_Before_All_Connected(i) = 500;
    end
    % Compute the gap statistic, using the 1-standard error method:
    system('python3 compute_gap_statistic.py'); % CODE VALID FOR LHS, EFAST, EXPERIMENTS 1-159, 255 -
%     system('python3 compute_gap_statistic_fixed.py'); % CODE VALID FOR EXPERIMENTS 160 - 225...
    Summary_Statistics.Gap_Statistic(i) = csvread('gap_statistic_puncta_model.txt');
    Summary_Statistics.Gap_Statistic_Standard_1_Error(i) = csvread('gap_std_1_error_statistic_puncta_model.txt');
    fprintf("\n\n\n");
end
% All infinite PH radii are >= 100 microns (largest filtration in ripser
% code)
% Summary_Statistics.Persistence_Homology_Largest_Radius_Before_All_Connected(isinf(Summary_Statistics.Persistence_Homology_Largest_Radius_Before_All_Connected)) = 100; % CODE VALID FOR LHS, EFAST, EXPERIMENTS 1-159
Summary_Statistics.Persistence_Homology_Largest_Radius_Before_All_Connected(isinf(Summary_Statistics.Persistence_Homology_Largest_Radius_Before_All_Connected)) = 500; % CODE VALID FOR LHS, EFAST, EXPERIMENTS 160 - ...

%% Record all average/standard deviation data over multiple realizations:
% Global = Average of cells in ALL realizations (not average of averages)
% Get global order parameters over all realizations:
global_average_order_parameter_all_cells = vecnorm(sum(Sum_of_All_Unit_Velocities_for_Order_Parameter_Calculation))./sum(Summary_Statistics.Number_of_Leaders + Summary_Statistics.Number_of_Followers);
global_average_order_parameter_leaders = vecnorm(sum(Sum_of_Leader_Unit_Velocities_for_Order_Parameter_Calculation))./sum(Summary_Statistics.Number_of_Leaders);
if sum(Summary_Statistics.Number_of_Followers) > 0
    global_average_order_parameter_followers = vecnorm(sum(Sum_of_Follower_Unit_Velocities_for_Order_Parameter_Calculation))./sum(Summary_Statistics.Number_of_Followers);
else
    global_average_order_parameter_followers = 0;
end
% Get global FN orientation over all realizations:
if sum(Summary_Statistics.Number_of_Oriented_FN_Puncta) > 0
    global_average_FN_orientation_angle = sum(Summary_Statistics.Average_FN_Puncta_Orientations)./sum(Summary_Statistics.Number_of_Oriented_FN_Puncta);
else
    global_average_FN_orientation_angle = 0;
end
Summary_Statistics.Average_FN_Puncta_Orientations = Summary_Statistics.Average_FN_Puncta_Orientations./Summary_Statistics.Number_of_Oriented_FN_Puncta;
Summary_Statistics.Average_FN_Puncta_Orientations(isnan(Summary_Statistics.Average_FN_Puncta_Orientations)) = 0;
Summary_Statistics.Average_FN_Puncta_Orientations(isinf(Summary_Statistics.Average_FN_Puncta_Orientations)) = 0;

% Get global statistic on average speed:
global_ave_speed = sum(Summary_Statistics.Average_Speed)./sum(Summary_Statistics.Number_of_Leaders + Summary_Statistics.Number_of_Followers);
% global_std_dev_speed = sqrt(1./(sum(Summary_Statistics.Number_of_Leaders + Summary_Statistics.Number_of_Followers)-1).*sum( sum(Speed_Second_Moment)-2*global_ave_speed*sum(Summary_Statistics.Average_Speed) + global_ave_speed.^2.*length(Speed_Second_Moment) ));
Summary_Statistics.Average_Speed = Summary_Statistics.Average_Speed./(Summary_Statistics.Number_of_Leaders + Summary_Statistics.Number_of_Followers);
Summary_Statistics.Global_Average_Leader_Speed = Summary_Statistics.Global_Average_Leader_Speed./sum(Summary_Statistics.Number_of_Leaders + Summary_Statistics.Number_of_Followers);
Summary_Statistics.Global_Average_Follower_Speed = Summary_Statistics.Global_Average_Follower_Speed./sum(Summary_Statistics.Number_of_Leaders + Summary_Statistics.Number_of_Followers);

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
% Create a csv file with all of the summary statistics:
csv_header = {'ID', 'average_estimated_stream_length (um)', 'standard_deviation_estimated_stream_length (um)',...%
    'average_estimated_stream_width (um)', 'standard_deviation_estimated_stream_width (um)',...%
    'average_max_stream_x_position (um)', 'standard_deviation_max_stream_x_position (um)', ... %
    'average_max_stream_y_width (um)', 'standard_deviation_max_stream_y_width (um)',...%
    'global_average_FN_puncta_orientation_over_all_realizations (rad)', 'average_FN_puncta_orientation_per_realization (rad)', 'standard_deviation_average_puncta_orientation_per_realization (rad)',...%
    'global_average_cell_speed_over_all_realizations (um/min)',...%
    'average_cell_speed_per_realization (um/min)','standard_deviation_cell_speed_per_realization (um/min)',...%
    'global_average_distance_to_nearest_cell_over_all_realizations (um)',...%
    'average_distance_to_nearest_cell_per_realization (um)', 'standard_deviation_average_distance_to_nearest_cell_per_realization (um)',...%
    'global_average_number_of_cells_in_50_micron_ball_across_all_realizations',...%
    'average_number_of_cells_in_50_micron_ball_per_realization','standard_deviation_average_number_of_cells_in_50_micron_ball_per_realization',...%
    'average_persistent_homology_radius_until_connected_stream (um)','standard_deviation_persistent_homology_radius_until_connected_stream (um)',...%
    'average_gap_statistic_value_over_all_realizations','standard_deviation_gap_statistic_value_over_all_realizations',...%
    'average_gap_statistic_standard_1_error_over_all_realizations','standard_deviation_gap_statistic_standard_1_error_over_all_realizations',...%
    'global_average_leader_x_position_at_final_time (um)','mean_average_leader_x_position_at_final_time', 'standard_deviation_average_leader_x_position_at_final_time'...%
    'global_average_follower_x_position_at_final_time (um)','mean_average_follower_x_position_at_final_time', 'standard_deviation_average_follower_x_position_at_final_time',...%
    'global_average_cell_order_parameter_over_all_realizations','average_cell_order_parameter_per_realization','standard_deviation_cell_order_parameter_per_realization',...%
    'global_average_leader_order_parameter_over_all_realizations','average_leader_order_parameter_per_realization','standard_deviation_leader_order_parameter_per_realization',...%
    'global_average_follower_order_parameter_over_all_realizations','average_follower_order_parameter_per_realization','standard_deviation_follower_order_parameter_per_realization',...%
    'average_number_of_leaders', 'standard_deviation_number_of_leaders',...%
    'average_number_of_followers', 'standard_deviation_number_of_followers',...%
    'average_number_of_oriented_FN_puncta', 'standard_deviation_number_of_oriented_FN_puncta'};%

%     'average_leader_angle (rad)', 'average_standard_deviation_leader_angle (rad)',...
%     'average_follower_angle (rad)', 'average_standard_deviation_follower_angle (rad)',...

    csv_header = strjoin(csv_header, ',');
%     temp_at_radius_50 = sum(Summary_Statistics.Histogram_Average_Cells_within_Radius_Statistic.*sum(cells_in_every_bin_per_realization,3),1)./sum(cells_in_every_bin_per_realization,3);
%     temp_at_radius_50 = temp_at_radius_50(10);
    % You already took into account the sum in the numerator for the
    % following fraction:
    mu_leaders = sum(Summary_Statistics.Leader_Cell_Average_x_position)./sum(Summary_Statistics.Number_of_Leaders);
    Summary_Statistics.Leader_Cell_Average_x_position = Summary_Statistics.Leader_Cell_Average_x_position./Summary_Statistics.Number_of_Leaders;
%     sigma_leaders = sqrt(1./sum(Summary_Statistics.Number_of_Leaders-1).*sum( sum(Leader_x_position_Second_Moment)- 2*mu_leaders*sum(Summary_Statistics.Leader_Cell_Average_x_position.*Summary_Statistics.Number_of_Leaders) + mu_leaders.^2.*length(Summary_Statistics.Number_of_Leaders)));
    if sum(Summary_Statistics.Number_of_Followers)>0
        % You already took into account the sum in the numerator for the
        % following fraction:
        mu_followers = sum(Summary_Statistics.Follower_Cell_Average_x_position)./sum(Summary_Statistics.Number_of_Followers);
%         sigma_followers = sqrt(1./sum(Summary_Statistics.Number_of_Followers-1).*sum( sum(Follower_x_position_Second_Moment)- 2*mu_leaders*sum(Summary_Statistics.Follower_Cell_Average_x_position.*Summary_Statistics.Number_of_Followers) + mu_followers.^2.*length(Summary_Statistics.Number_of_Followers)));
        Summary_Statistics.Follower_Cell_Average_x_position = Summary_Statistics.Follower_Cell_Average_x_position./Summary_Statistics.Number_of_Followers;
    else
        mu_followers = 0;
%         sigma_followers = 0;
    end
    csv_values = [parameter_set_number, mean(Summary_Statistics.Est_Stream_Length), std(Summary_Statistics.Est_Stream_Length),...
        mean(Summary_Statistics.Est_Stream_Width), std(Summary_Statistics.Est_Stream_Width),...
        mean(Summary_Statistics.Max_X_Length_of_Stream), std(Summary_Statistics.Max_X_Length_of_Stream),...
        mean(Summary_Statistics.Max_Y_Length_of_Stream), std(Summary_Statistics.Max_Y_Length_of_Stream),...
        global_average_FN_orientation_angle, mean(Summary_Statistics.Average_FN_Puncta_Orientations), std(Summary_Statistics.Average_FN_Puncta_Orientations),... % first is true mean, second is mean and std dev of averages
        global_ave_speed,... % true mean
        mean(Summary_Statistics.Average_Speed),std(Summary_Statistics.Average_Speed),... % mean and stdev of averages
        sum(Summary_Statistics.Average_Min_Distance_per_Realization.*(Summary_Statistics.Number_of_Leaders + Summary_Statistics.Number_of_Followers))./sum(Summary_Statistics.Number_of_Leaders + Summary_Statistics.Number_of_Followers),... %true mean
        mean(Summary_Statistics.Average_Min_Distance_per_Realization), std(Summary_Statistics.Average_Min_Distance_per_Realization),...% mean and stdev of averages
        sum(Summary_Statistics.Average_Number_of_Cells_within_50_um_Ball_per_Realization.*(Summary_Statistics.Number_of_Leaders + Summary_Statistics.Number_of_Followers))./sum(Summary_Statistics.Number_of_Leaders + Summary_Statistics.Number_of_Followers),... %true mean
        mean(Summary_Statistics.Average_Number_of_Cells_within_50_um_Ball_per_Realization), std(Summary_Statistics.Average_Number_of_Cells_within_50_um_Ball_per_Realization),... % mean and stdev of averages
        mean(Summary_Statistics.Persistence_Homology_Largest_Radius_Before_All_Connected), std(Summary_Statistics.Persistence_Homology_Largest_Radius_Before_All_Connected),...
        mean(Summary_Statistics.Gap_Statistic), std(Summary_Statistics.Gap_Statistic),...
        mean(Summary_Statistics.Gap_Statistic_Standard_1_Error), std(Summary_Statistics.Gap_Statistic_Standard_1_Error),...
        mu_leaders, mean(Summary_Statistics.Leader_Cell_Average_x_position), std(Summary_Statistics.Leader_Cell_Average_x_position),...
        mu_followers, mean(Summary_Statistics.Follower_Cell_Average_x_position), std(Summary_Statistics.Follower_Cell_Average_x_position),...
        global_average_order_parameter_all_cells, mean(Summary_Statistics.Order_Parameter_All_Cells), std(Summary_Statistics.Order_Parameter_All_Cells),...
        global_average_order_parameter_leaders, mean(Summary_Statistics.Order_Parameter_Leaders), std(Summary_Statistics.Order_Parameter_Leaders),...
        global_average_order_parameter_followers, mean(Summary_Statistics.Order_Parameter_Followers),std(Summary_Statistics.Order_Parameter_Followers),...
        mean(Summary_Statistics.Number_of_Leaders), std(Summary_Statistics.Number_of_Leaders),...
        mean(Summary_Statistics.Number_of_Followers), std(Summary_Statistics.Number_of_Followers),...
        mean(Summary_Statistics.Number_of_Oriented_FN_Puncta), std(Summary_Statistics.Number_of_Oriented_FN_Puncta)];
        
%         mean(temp_at_radius_50),...
%         std(temp_at_radius_50),...

  
    
%         atan2(sum(sin(Summary_Statistics.Average_Leader_Angle)), sum(cos(Summary_Statistics.Average_Leader_Angle))),...
%         mean(Summary_Statistics.Standard_Deviation_Leader_Angle),...
%         atan2(sum(sin(Summary_Statistics.Average_Follower_Angle)), sum(cos(Summary_Statistics.Average_Follower_Angle))),...
%         mean(Summary_Statistics.Standard_Deviation_Leader_Angle),...

%% Save the config file in the output folder
save('./output/summary_statistics.mat', 'Summary_Statistics');
system('cp ./config/PhysiCell_settings.xml ./output');
%% Make a movie of the final realization. 
video_title = num2str(parameter_set_number, '%04.f');
folder_location = './output/';
make_video = true;
if make_video
    make_ncc_model_puncta_movie_experiments(video_title, folder_location);

   %% Transfer Everything to a New Directory and compress the resulting folder:
   system(strcat('mkdir ./NCC_puncta_model_output_data_experiments/', video_title)); 
   system(strcat('mv ./output/* ./NCC_puncta_model_output_data_experiments/',video_title) );
   system('touch ./output/empty.txt');
   system(['cd ./NCC_puncta_model_output_data_experiments/',video_title,' && tar -zcf ',video_title,'.tar.gz * && mv ',video_title,'.tar.gz ../ && cd -']);
   system(['rm -rf ./NCC_puncta_model_output_data_experiments/',video_title]);
end
end