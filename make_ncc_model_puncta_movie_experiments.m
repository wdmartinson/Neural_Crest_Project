function make_ncc_model_puncta_movie_experiments(video_title, folder_location)
str = computer;
if str(1) == 'M'
    setenv('PATH', '/usr/local/bin:/usr/local/sbin:/usr/bin:/bin:/usr/sbin:/sbin');
end
close all;
% Find out how many XML files there are
listing = dir(strcat(folder_location,'*.xml'));
N = length(listing);

% Load the first MCDS File
MCDS = read_MultiCellDS_xml(listing(2).name, folder_location);
dx = MCDS.mesh.X_coordinates(2)-MCDS.mesh.X_coordinates(1);
dy = MCDS.mesh.Y_coordinates(2)-MCDS.mesh.Y_coordinates(1);
Xmin = MCDS.mesh.X_coordinates(1) - dx/2;
Ymin = MCDS.mesh.Y_coordinates(1) - dy/2;
Xmax = MCDS.mesh.X_coordinates(end) + dx/2;
Ymax = MCDS.mesh.Y_coordinates(end) + dy/2;

MCDS = read_MultiCellDS_xml(listing(end).name, folder_location);
FN_max = max(max(MCDS.continuum_variables(1).data));

fig1 = figure('units', 'inches', 'position', [0,0,7,5]);
    imagesc('xdata', MCDS.mesh.X_coordinates, 'ydata', fliplr(MCDS.mesh.Y_coordinates), 'cdata', NaN);
    colormap([1 1 1; 0 0 1]);
    xlabel('$x$, $\mu$m', 'interpreter', 'latex', 'fontsize', 16);
    ylabel('$y$, $\mu$m', 'interpreter', 'latex', 'fontsize', 16);
    
    axis square;
    set(gca, 'nextplot', 'replacechildren', 'visible','on');
    ax1 = gca;
    title('','interpreter', 'latex', 'fontsize', 20);
    set(get(ax1, 'title'), 'string', ['$t$ = ',num2str(MCDS.metadata.current_time), ' ', MCDS.metadata.time_units]);
%     cbar = colorbar;
    caxis([0 FN_max]);
%     set(get(cbar, 'title'), 'string', 'Fibronectin, $\mu$M', 'interpreter', 'latex', 'fontsize', 16);
    ax1.Title.Interpreter = 'latex';
    ax1.FontSize = 16;
    xlim([Xmin, Xmax]);
    ylim([Ymin, Ymax]);
    
    
fig2 = figure('units', 'inches', 'position', [0,0,7,5]);
    xlabel('$x$, $\mu$m', 'interpreter', 'latex', 'fontsize', 16);
    ylabel('$y$, $\mu$m', 'interpreter', 'latex', 'fontsize', 16);
    
    axis square; box on;
    set(gca, 'nextplot', 'replacechildren', 'visible','on');
    ax2 = gca;
    title('','interpreter', 'latex', 'fontsize', 20);
    set(get(ax2, 'title'), 'string', ['$t$ = ',num2str(MCDS.metadata.current_time), ' ', MCDS.metadata.time_units]);
    ax2.Title.Interpreter = 'latex';
    ax2.FontSize = 16;
    xlim([Xmin, Xmax]);
    ylim([Ymin, Ymax]);

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
    
a_name = video_title;
b_name = [video_title,'-FN_orientation'];
c_name = [video_title,'-cells_with_directional_cues'];
a = VideoWriter(a_name);
a.Quality = 100;
a.FrameRate = 5;
b = VideoWriter(b_name);
b.Quality = 100;
b.FrameRate = 5;
c = VideoWriter(c_name);
c.Quality = 100;
c.FrameRate = 5;
open(a);
open(b);
open(c);
% Load in the MCDS files. The first three files are 'config.xml', 'final.xml', and
% 'initial.xml', so start from i = 4:
for i = 4:N
    figure(fig1);
    cla;
    % Load the MCDS file:
    MCDS = read_MultiCellDS_xml(listing(i).name, folder_location);
    % Extract the cells
    cells = MCDS.discrete_cells;
    Leaders = (cells.metadata.type == 0);
    Followers = (cells.metadata.type == 1);
    Puncta = (cells.metadata.type == 2);
    % Extract the microenvironment:
    FN_field = MCDS.continuum_variables(1).data;
    % Plot the FN field and then the cells themselves as a scatter
    % plot:
    leader_radii = (cells.phenotype.geometrical_properties.volumes.nuclear(Leaders)./4./pi.*3).^(1/3);
    follower_radii = (cells.phenotype.geometrical_properties.volumes.nuclear(Followers)./4./pi.*3).^(1/3);
    imagesc('xdata', MCDS.mesh.X_coordinates, 'ydata', fliplr(MCDS.mesh.Y_coordinates), 'cdata', flipud(FN_field));
%     cbar = colorbar; 
    caxis([0 FN_max]);
    hold on;
    h1 = viscircles([cells.state.position(Leaders,1) cells.state.position(Leaders,2)], leader_radii, 'color', 'k', 'linewidth', 0.75);
    if ~isempty(cells.state.position(Leaders,1))
        x_leaders = h1.Children(1).XData;
        y_leaders = h1.Children(1).YData;
        idx = [0, find(isnan(x_leaders))];
        for ii = 1:length(idx)-1
            fill(x_leaders(idx(ii)+1:idx(ii+1)-1), y_leaders(idx(ii)+1:idx(ii+1)-1), 'k', 'facealpha', 0.5);
        end
    end
    h2 = viscircles([cells.state.position(Followers,1) cells.state.position(Followers,2)], follower_radii, 'color', 'r', 'linewidth', 0.75);
    if ~isempty(cells.state.position(Followers,1))
        x_followers = h2.Children(1).XData;
        y_followers = h2.Children(1).YData;
        idx = [0, find(isnan(x_followers))];
        for ii = 1:length(idx)-1
            fill(x_followers(idx(ii)+1:idx(ii+1)-1), y_followers(idx(ii)+1:idx(ii+1)-1), 'r', 'facealpha', 0.5);
        end
    end
    Puncta1 = logical(Puncta.*cells.custom.time_since_last_filopodia_drop);
    quiver(cells.state.position(Leaders, 1), cells.state.position(Leaders, 2), ...
        15*cells.custom.total_velocity(Leaders, 1)./sqrt(sum(cells.custom.total_velocity(Leaders, :).^2, 2)), 15*cells.custom.total_velocity(Leaders,2)./sqrt(sum(cells.custom.total_velocity(Leaders, :).^2, 2)),...
        0, 'color', 'k');
    quiver(cells.state.position(Followers, 1), cells.state.position(Followers, 2), ...
        15*cells.custom.total_velocity(Followers, 1)./sqrt(sum(cells.custom.total_velocity(Followers, :).^2, 2)), 15*cells.custom.total_velocity(Followers,2)./sqrt(sum(cells.custom.total_velocity(Followers, :).^2, 2)),...
        0, 'color', 'r');
     quiver(cells.state.position(Puncta1, 1), cells.state.position(Puncta1, 2), ...
        10*cos(cells.custom.time_until_next_filopodia_drop(Puncta1)'), 10*sin(cells.custom.time_until_next_filopodia_drop(Puncta1)'),...
        0, 'color', 'b');
    % Housekeeping: Keeping the video quality the same throughout
    xlim([Xmin, Xmax]); ylim([Ymin, Ymax]);
    set(get(ax1, 'title'), 'string', ['$t$ = ',num2str(MCDS.metadata.current_time), ' ', MCDS.metadata.time_units]);
%     set(get(cbar, 'title'), 'string', 'Fibronectin, $\mu$M', 'interpreter', 'latex', 'fontsize', 16);
    ax1.Title.Interpreter = 'latex';
    ax1.FontSize = 20;
    xlabel('$x$, $\mu$m', 'interpreter', 'latex', 'fontsize', 16);
    ylabel('$y$, $\mu$m', 'interpreter', 'latex', 'fontsize', 16); 
    hold off;
    
    figure(fig2);
    cla;
    hold on;
    quiver(cells.state.position(Puncta1, 1), cells.state.position(Puncta1, 2), ...
        10*cos(cells.custom.time_until_next_filopodia_drop(Puncta1)'), 10*sin(cells.custom.time_until_next_filopodia_drop(Puncta1)'),...
        0, 'color', 'b');
    xlim([Xmin, Xmax]); ylim([Ymin, Ymax]);
    set(get(ax2, 'title'), 'string', ['$t$ = ',num2str(MCDS.metadata.current_time), ' ', MCDS.metadata.time_units]);
    ax2.Title.Interpreter = 'latex';
    ax2.FontSize = 20;
    xlabel('$x$, $\mu$m', 'interpreter', 'latex', 'fontsize', 16);
    ylabel('$y$, $\mu$m', 'interpreter', 'latex', 'fontsize', 16); 
    hold off;
    
    figure(fig3);
    cla;
    hold on;
    h1 = viscircles([cells.state.position(Leaders,1) cells.state.position(Leaders,2)], leader_radii, 'color', 'k', 'linewidth', 0.75);
    if ~isempty(cells.state.position(Leaders,1))
        x_leaders = h1.Children(1).XData;
        y_leaders = h1.Children(1).YData;
        idx = [0, find(isnan(x_leaders))];
        for ii = 1:length(idx)-1
            fill(x_leaders(idx(ii)+1:idx(ii+1)-1), y_leaders(idx(ii)+1:idx(ii+1)-1), 'k', 'facealpha', 0.5);
        end
    end
    h2 = viscircles([cells.state.position(Followers,1) cells.state.position(Followers,2)], follower_radii, 'color', 'r', 'linewidth', 0.75);
    if ~isempty(cells.state.position(Followers,1))
        x_followers = h2.Children(1).XData;
        y_followers = h2.Children(1).YData;
        idx = [0, find(isnan(x_followers))];
        for ii = 1:length(idx)-1
            fill(x_followers(idx(ii)+1:idx(ii+1)-1), y_followers(idx(ii)+1:idx(ii+1)-1), 'r', 'facealpha', 0.5);
        end
    end
    quiver(cells.state.position(Leaders, 1), cells.state.position(Leaders, 2), ...
        20*cells.custom.total_velocity(Leaders, 1)./sqrt(sum(cells.custom.total_velocity(Leaders, :).^2, 2)), 20*cells.custom.total_velocity(Leaders,2)./sqrt(sum(cells.custom.total_velocity(Leaders, :).^2, 2)),...
        0, 'color', 'k');
    quiver(cells.state.position(Followers, 1), cells.state.position(Followers, 2), ...
        20*cells.custom.total_velocity(Followers, 1)./sqrt(sum(cells.custom.total_velocity(Followers, :).^2, 2)), 20*cells.custom.total_velocity(Followers,2)./sqrt(sum(cells.custom.total_velocity(Followers, :).^2, 2)),...
        0, 'color', 'k');
    quiver(cells.state.position(Leaders, 1), cells.state.position(Leaders, 2), ...
        10*cells.custom.haptotaxis_direction(Leaders, 1), 10*cells.custom.haptotaxis_direction(Leaders,2),...
        0, 'color', 'r');
    quiver(cells.state.position(Followers, 1), cells.state.position(Followers, 2), ...
        10*cells.custom.haptotaxis_direction(Followers, 1), 10*cells.custom.haptotaxis_direction(Followers,2),...
        0, 'color', 'r');
    quiver(cells.state.position(Leaders, 1), cells.state.position(Leaders, 2), ...
        10*cells.custom.contact_guidance_direction(Leaders, 1), 10*cells.custom.contact_guidance_direction(Leaders,2),...
        0, 'color', [.5 0 .5]);
    quiver(cells.state.position(Followers, 1), cells.state.position(Followers, 2), ...
        10*cells.custom.contact_guidance_direction(Followers, 1), 10*cells.custom.contact_guidance_direction(Followers,2),...
        0, 'color', [.5 0 .5]);
    if ~isfield(cells.custom, 'b_adhesion')
        quiver(cells.state.position(Leaders, 1), cells.state.position(Leaders, 2), ...
            10*cells.custom.volume_exclusion_direction(Leaders, 1), 10*cells.custom.volume_exclusion_direction(Leaders,2),...
            0, 'color', 'b');
        quiver(cells.state.position(Followers, 1), cells.state.position(Followers, 2), ...
            10*cells.custom.volume_exclusion_direction(Followers, 1), 10*cells.custom.volume_exclusion_direction(Followers,2),...
            0, 'color', 'b');
    else
        quiver(cells.state.position(Leaders, 1), cells.state.position(Leaders, 2), ...
            10.*sign(cells.custom.b_vol_exclusion(Leaders)-cells.custom.b_adhesion(Leaders))'.*cells.custom.volume_exclusion_direction(Leaders, 1), 10.*sign(cells.custom.b_vol_exclusion(Leaders)-cells.custom.b_adhesion(Leaders))'.*cells.custom.volume_exclusion_direction(Leaders,2),...
            0, 'color', 'b');
        quiver(cells.state.position(Followers, 1), cells.state.position(Followers, 2), ...
            10.*sign(cells.custom.b_vol_exclusion(Followers)-cells.custom.b_adhesion(Followers))'.*cells.custom.volume_exclusion_direction(Followers, 1), 10.*sign(cells.custom.b_vol_exclusion(Followers)-cells.custom.b_adhesion(Followers))'.*cells.custom.volume_exclusion_direction(Followers,2),...
            0, 'color', 'b');
    end
    xlim([Xmin, Xmax]); ylim([Ymin, Ymax]);
    set(get(ax3, 'title'), 'string', ['$t$ = ',num2str(MCDS.metadata.current_time), ' ', MCDS.metadata.time_units]);
    ax3.Title.Interpreter = 'latex';
    ax3.FontSize = 20;
    xlabel('$x$, $\mu$m', 'interpreter', 'latex', 'fontsize', 16);
    ylabel('$y$, $\mu$m', 'interpreter', 'latex', 'fontsize', 16); 
    hold off;
    
    % Take the snapshot
    frame = getframe(fig1);
    frame2 = getframe(fig2);
    frame3 = getframe(fig3);
    writeVideo(a, frame);
    writeVideo(b, frame2);
    writeVideo(c, frame3);
end % for i
close(a);
close(b);
close(c);

% Convert the AVI file to an MP4:
system(['ffmpeg -n -i ',a_name,'.avi ',a_name,'.mp4']);
system(['ffmpeg -n -i ',b_name,'.avi ',b_name,'.mp4']);
system(['ffmpeg -n -i ',c_name,'.avi ',c_name,'.mp4']);

% Delete the AVI file:
system(['rm -f ', a_name,'.avi']);
system(['rm -f ', b_name,'.avi']);
system(['rm -f ', c_name,'.avi']);

% Move the MP4 to the storage folder:
system(['mv ',a_name,'.mp4 ',folder_location,a_name,'.mp4']);
system(['mv ',b_name,'.mp4 ',folder_location,b_name,'.mp4']);
system(['mv ',c_name,'.mp4 ',folder_location,c_name,'.mp4']);

leaders_and_followers = logical(Leaders+Followers);

% Simple scatter plot at final time point, with unit vector velocities, for zegami:
saveas(fig1,[folder_location,video_title,'.svg'],'svg');
saveas(fig1,[folder_location,video_title,'.png'], 'png');

% Also plot just the fibronectin field, for posterity's sake:
saveas(fig2,[folder_location,video_title,'-FN_orientation.svg'], 'svg');
saveas(fig2,[folder_location,video_title,'-FN_orientation.png'], 'png');

% And save the last directional cue plot, for posterity's sake:
saveas(fig3,[folder_location,video_title,'-cells_with_directional_cues.svg'], 'svg');
saveas(fig3,[folder_location,video_title,'-cells_with_directional_cues.png'], 'png');

load([folder_location, 'summary_statistics.mat'], 'Summary_Statistics');

if isfield(Summary_Statistics, 'Histogram_Bin_Edges')
   midpoints = 0.5*(Summary_Statistics.Histogram_Bin_Edges(1:end-1) + Summary_Statistics.Histogram_Bin_Edges(2:end)); 
   figure('units', 'inches', 'position', [0,0,7,5]);
   p1 = bar(midpoints, Summary_Statistics.Histogram_Average_Leader_Cell_x_Position, 'FaceAlpha', 0.5);
   hold on;
   p2 = bar(midpoints, Summary_Statistics.Histogram_Average_Follower_Cell_x_Position, 'FaceAlpha', 0.5);
   legend([p1(1) p2(1)], 'Leader Cell', 'Follower Cell');
   xlabel('$x$, $\mu$m', 'interpreter', 'latex', 'fontsize', 16);
   ylabel('Average Number of Cells', 'interpreter', 'latex', 'fontsize', 16);
   title('Average Cell Distribution in $x$-direction', 'interpreter', 'latex', 'fontsize', 16);
   saveas(gcf, [folder_location,'a-Final_Histogram_x_position.svg'], 'svg');
end

if isfield(Summary_Statistics, 'Est_Stream_Length')
    figure('units', 'inches', 'position', [0,0,7,5]);
    p1 = plot(1:length(Summary_Statistics.Est_Stream_Length), Summary_Statistics.Est_Stream_Length, 'ob');
    hold on;
    p2 = fplot(@(x) 0.*x + mean(Summary_Statistics.Est_Stream_Length), [1, length(Summary_Statistics.Est_Stream_Length)], '--k');
    legend([p1(1), p2(1)], 'Radius for Particular Realization', 'Average Radius for all Realizations');
    xlabel('Realization', 'interpreter', 'latex', 'fontsize', 16);
    ylabel('Length, $\mu$m', 'interpreter', 'latex', 'fontsize', 16);
    title('Estimated Stream Length', 'interpreter', 'latex', 'fontsize', 16);
    saveas(gcf, [folder_location,'b-Estimated_Stream_Length_in_FN_orientation.svg'], 'svg');
end

if isfield(Summary_Statistics, 'Est_Stream_Width')
   figure('units', 'inches', 'position', [0,0,7,5]);
   p1 = plot((1:length(Summary_Statistics.Est_Stream_Width))', Summary_Statistics.Est_Stream_Width, 'ob');
   hold on;
   p2 = fplot(@(x) 0*x + mean(Summary_Statistics.Est_Stream_Width), [1, length(Summary_Statistics.Est_Stream_Width)], '--k');
   legend([p1(1), p2(1)], 'Statistic for Particular Realization', 'Mean for all Realizations');
   xlabel('Realization', 'interpreter', 'latex', 'fontsize', 16);
   ylabel('Width, $\mu m$', 'interpreter', 'latex', 'fontsize', 16);
   title('Estimated Stream Width', 'interpreter', 'latex', 'fontsize', 16);
   saveas(gcf, [folder_location,'c-Estimated_Stream_Width_in_FN_orientation.svg'], 'svg');
end

if isfield(Summary_Statistics, 'Average_FN_Puncta_Orientations')
   figure('units', 'inches', 'position', [0,0,7,5]);
%    p1 = plot((1:length(Summary_Statistics.Average_FN_Puncta_Orientations))', Summary_Statistics.Average_FN_Puncta_Orientations, 'ob');
%    hold on;
%    p2 = fplot(@(x) 0*x + mean(Summary_Statistics.Average_FN_Puncta_Orientations), [1, length(Summary_Statistics.Average_FN_Puncta_Orientations)], '--k');
%    ylim([-pi, pi]);
%    xlabel('Realization', 'interpreter', 'latex', 'fontsize', 16);
%    ylabel('Angle, (rad)', 'interpreter', 'latex', 'fontsize', 16);
%    legend([p1(1), p2(1)], 'Statistic for Particular Realization', 'Mean for all Realizations');
   p1 = polarscatter(Summary_Statistics.Average_FN_Puncta_Orientations, (1:length(Summary_Statistics.Average_FN_Puncta_Orientations))', 'ob');
   hold on;
   p2 = polarplot(mean(Summary_Statistics.Average_FN_Puncta_Orientations).*ones(length(Summary_Statistics.Average_FN_Puncta_Orientations)+1, 1), (0:length(Summary_Statistics.Average_FN_Puncta_Orientations))', '--k');
   legend([p1(1), p2(1)], 'Statistic for Particular Realization', 'Mean for all Realizations');
   title('Average FN Puncta Orientation', 'interpreter', 'latex', 'fontsize', 16);
   saveas(gcf, [folder_location,'d-Average_FN_Puncta_Orientation.svg'], 'svg');
end

if isfield(Summary_Statistics, 'Max_X_Length_of_Stream')
   figure('units', 'inches', 'position', [0,0,7,5]);
   p1 = plot((1:length(Summary_Statistics.Max_X_Length_of_Stream))', Summary_Statistics.Max_X_Length_of_Stream, 'ob');
   hold on;
   p2 = fplot(@(x) 0*x + mean(Summary_Statistics.Max_X_Length_of_Stream), [1, length(Summary_Statistics.Max_X_Length_of_Stream)], '--k');
   legend([p1(1), p2(1)], 'Statistic for Particular Realization', 'Mean for all Realizations');
   xlabel('Realization', 'interpreter', 'latex', 'fontsize', 16);
   ylabel('Distance, $\mu m$', 'interpreter', 'latex', 'fontsize', 16);
   title('Max X-Length of Stream', 'interpreter', 'latex', 'fontsize', 16);
   saveas(gcf, [folder_location,'e-Max_X_Length_of_Stream.svg'], 'svg');
end

if isfield(Summary_Statistics, 'Max_Y_Length_of_Stream')
   figure('units', 'inches', 'position', [0,0,7,5]);
   p1 = plot((1:length(Summary_Statistics.Max_Y_Length_of_Stream))', Summary_Statistics.Max_Y_Length_of_Stream, 'ob');
   hold on;
   p2 = fplot(@(x) 0*x + mean(Summary_Statistics.Max_Y_Length_of_Stream), [1, length(Summary_Statistics.Max_Y_Length_of_Stream)], '--k');
   legend([p1(1), p2(1)], 'Statistic for Particular Realization', 'Mean for all Realizations');
   xlabel('Realization', 'interpreter', 'latex', 'fontsize', 16);
   ylabel('Distance, $\mu m$', 'interpreter', 'latex', 'fontsize', 16);
   title('Max Y-Width of Stream', 'interpreter', 'latex', 'fontsize', 16);
   saveas(gcf, [folder_location,'f-Max_Y_Width_of_Stream.svg'], 'svg');
end

if isfield(Summary_Statistics, 'Average_Min_Distance_per_Realization')
   figure('units', 'inches', 'position', [0,0,7,5]);
   p1 = plot((1:length(Summary_Statistics.Average_Min_Distance_per_Realization))', Summary_Statistics.Average_Min_Distance_per_Realization, 'ob');
   hold on;
   p2 = fplot(@(x) 0*x + mean(Summary_Statistics.Average_Min_Distance_per_Realization), [1, length(Summary_Statistics.Average_Min_Distance_per_Realization)], '--k');
   legend([p1(1), p2(1)], 'Statistic for Particular Realization', 'Mean for all Realizations');
   xlabel('Realization', 'interpreter', 'latex', 'fontsize', 16);
   ylabel('Distance, $\mu m$', 'interpreter', 'latex', 'fontsize', 16);
   title('Average Distance to Nearest Cell per Realization', 'interpreter', 'latex', 'fontsize', 16);
   saveas(gcf, [folder_location,'g-Average_Min_Distance_per_Realization.svg'], 'svg');
end

if isfield(Summary_Statistics, 'Histogram_Average_Minimum_Distance_Measure')
   midpoints = 0.5*(Summary_Statistics.Histogram_Bin_Edges(1:end-1) + Summary_Statistics.Histogram_Bin_Edges(2:end)); 
   figure('units', 'inches', 'position', [0,0,7,5]);
   bar(midpoints, Summary_Statistics.Histogram_Average_Minimum_Distance_Measure);
   xlabel('$x$, $\mu$m', 'interpreter', 'latex', 'fontsize', 16);
   ylabel('Distance, $\mu$m', 'interpreter', 'latex', 'fontsize', 16);
   title('Average distance to the nearest cell', 'interpreter', 'latex', 'fontsize', 16);
   saveas(gcf, [folder_location,'h-Histogram_Average_Minimum_Distance_Measure.svg'], 'svg');
end

if isfield(Summary_Statistics, 'Histogram_Average_Cells_within_Radius_Statistic')
    % Heat map of the persistence homology statistic
   midpoints = 0.5*(Summary_Statistics.Histogram_Bin_Edges(1:end-1) + Summary_Statistics.Histogram_Bin_Edges(2:end));
   figure('units', 'inches', 'position', [0,0,7,5]);
   imagesc('xdata', midpoints, 'ydata', Summary_Statistics.Density_Radii, 'cdata', Summary_Statistics.Histogram_Average_Cells_within_Radius_Statistic');
   axis tight;
   cbar = colorbar;
   xlabel('$x$, $\mu$m', 'interpreter', 'latex', 'fontsize', 16);
   ylabel('Radius $\epsilon$, $\mu$m', 'interpreter', 'latex', 'fontsize', 16);
   title('Number of cells within ball of radius $\epsilon$', 'interpreter', 'latex', 'fontsize', 16);
   set(get(cbar, 'title'), 'string', 'No. of Cells', 'interpreter', 'latex', 'fontsize', 16);
   saveas(gcf, [folder_location,'i-Histogram_Average_Cells_within_Radius_Heat_Map.svg'], 'svg');
   
   % Plot the results as lines for all values of r
   figure('units', 'inches', 'position', [0,0,7,5]);
   p1 = plot(midpoints, Summary_Statistics.Histogram_Average_Cells_within_Radius_Statistic);
   xlabel('$x$, $\mu$m', 'interpreter', 'latex', 'fontsize', 16);
   ylabel('Number of Cells', 'interpreter', 'latex', 'fontsize', 16);
   strings = cell(1, length(Summary_Statistics.Density_Radii));
   for i = 1:length(Summary_Statistics.Density_Radii)
    strings{i} = ['$\epsilon = ',num2str(Summary_Statistics.Density_Radii(i)),' \mu m$'];
   end
   legend(p1, strings, 'interpreter', 'latex');
   title('Number of cells within ball of radius $\epsilon$', 'interpreter', 'latex', 'fontsize', 16);
   saveas(gcf, [folder_location,'j-Histogram_Average_Cells_within_Radius_Line_Plot.svg'], 'svg');
   
   % Also plot the results for a particular value of r. Choose the halfway
   % one for simplicity
   figure('units', 'inches', 'position', [0,0,7,5]);
   bar(midpoints, Summary_Statistics.Histogram_Average_Cells_within_Radius_Statistic(:, 10), 'k');
   xlim([Xmin, Xmax]);
   xlabel('$x$, $\mu$m', 'interpreter', 'latex', 'fontsize', 16);
   ylabel('Number of Cells', 'interpreter', 'latex', 'fontsize', 16);
   title(['Radius = ', num2str(Summary_Statistics.Density_Radii(10)),' $\mu m$'] ,'interpreter', 'latex', 'fontsize', 16);
   saveas(gcf, [folder_location,'k-Histogram_Average_Cells_within_particular_radius.svg'], 'svg');
end

if isfield(Summary_Statistics, 'Persistence_Homology_Largest_Radius_Before_All_Connected')
   % Plot the radius found by persistence homology of every realization. Also plot the average radius across all realizations:
   figure('units', 'inches', 'position', [0,0,7,5]);
   p1 = plot(1:length(Summary_Statistics.Persistence_Homology_Largest_Radius_Before_All_Connected), Summary_Statistics.Persistence_Homology_Largest_Radius_Before_All_Connected, 'o');
   hold on;
   p2 = fplot(@(x) 0.*x + mean(Summary_Statistics.Persistence_Homology_Largest_Radius_Before_All_Connected), [1, length(Summary_Statistics.Persistence_Homology_Largest_Radius_Before_All_Connected)], '--k');
   legend([p1(1), p2(1)], 'Radius for Particular Realization', 'Average Radius for all Realizations');
   xlabel('Realization', 'interpreter', 'latex', 'fontsize', 16);
   ylabel('Radius ($\mu$m)', 'interpreter', 'latex', 'fontsize', 16);
   title('Radii found using Persistent Homology');
   saveas(gcf, [folder_location,'l-Actual_Persistent_Homology_Radii.svg'], 'svg');
   % Also check to make sure that the radius is correct for the final
   % realization:
   figure('units', 'inches', 'position', [0,0,7,5]);
   axis square;
   hold on; box on;

    viscircles([cells.state.position(Leaders,1) cells.state.position(Leaders,2)], leader_radii, 'color', 'k', 'linewidth', 0.5);
    viscircles([cells.state.position(Followers,1) cells.state.position(Followers,2)], follower_radii, 'color', 'r', 'linewidth', 0.5);
    quiver(cells.state.position(Leaders, 1), cells.state.position(Leaders, 2), ...
        cells.custom.total_velocity(Leaders, 1)./sqrt(sum(cells.custom.total_velocity(Leaders, :).^2, 2)), cells.custom.total_velocity(Leaders,2)./sqrt(sum(cells.custom.total_velocity(Leaders, :).^2, 2)),...
        'k', 'AutoScaleFactor', 0.2);
    quiver(cells.state.position(Followers, 1), cells.state.position(Followers, 2), ...
        cells.custom.total_velocity(Followers, 1)./sqrt(sum(cells.custom.total_velocity(Followers, :).^2, 2)), cells.custom.total_velocity(Followers,2)./sqrt(sum(cells.custom.total_velocity(Followers, :).^2, 2)),...
        'r', 'AutoScaleFactor', 0.2);
    viscircles(cells.state.position(leaders_and_followers, 1:2), Summary_Statistics.Persistence_Homology_Largest_Radius_Before_All_Connected(end).*ones(size(cells.state.position(leaders_and_followers, :), 1), 1), 'Color','b');
    xlim([Xmin, Xmax]);
    ylim([Ymin, Ymax]);
    xlabel('$x$, $\mu$m', 'interpreter', 'latex', 'fontsize', 16);
    ylabel('$y$, $\mu$m', 'interpreter', 'latex', 'fontsize', 16);   
    title(['$t$ = ',num2str(MCDS.metadata.current_time), ' ', MCDS.metadata.time_units],'interpreter', 'latex', 'fontsize', 16);
    saveas(gcf, [folder_location,'m-Check_Actual_Persistent_Homology_Radii.svg'], 'svg');
end

if isfield(Summary_Statistics, 'Gap_Statistic')
   figure('units', 'inches', 'position', [0,0,7,5]);
   p1 = plot(1:length(Summary_Statistics.Gap_Statistic), Summary_Statistics.Gap_Statistic, 'o');
   hold on;
   p2 = fplot(@(x) 0*x + mean(Summary_Statistics.Gap_Statistic), [1, length(Summary_Statistics.Gap_Statistic)], '--k');
   legend([p1(1), p2(1)], 'Statistic for Particular Realization', 'Average Statistic for all Realizations');
   xlabel('Realization', 'interpreter', 'latex', 'fontsize', 16);
   ylabel('Number of Clusters', 'interpreter', 'latex', 'fontsize', 16);
   title('Gap Statistic');
   saveas(gcf, [folder_location,'n-Gap_Statistic.svg'], 'svg');
   
   figure('units', 'inches', 'position', [0,0,7,5]);
   p1 = plot(1:length(Summary_Statistics.Gap_Statistic_Standard_1_Error), Summary_Statistics.Gap_Statistic_Standard_1_Error, 'o');
   hold on;
   p2 = fplot(@(x) 0*x + mean(Summary_Statistics.Gap_Statistic_Standard_1_Error), [1, length(Summary_Statistics.Gap_Statistic_Standard_1_Error)], '--k');
   legend([p1(1), p2(1)], 'Statistic for Particular Realization', 'Average Statistic for all Realizations');
   xlabel('Realization', 'interpreter', 'latex', 'fontsize', 16);
   ylabel('Number of Clusters', 'interpreter', 'latex', 'fontsize', 16);
   title('Gap Statistic, Standard 1 Error', 'interpreter', 'latex', 'fontsize', 16);
   saveas(gcf, [folder_location,'o-Gap_Statistic_Standard_1_Error.svg'], 'svg');
   
   figure('units', 'inches', 'position', [0,0,7,5]);
    hold on;
    pos = cells.state.position(leaders_and_followers,:);
    radii = (cells.phenotype.geometrical_properties.volumes.nuclear(leaders_and_followers)./4./pi.*3).^(1/3);
    idx = kmeans(pos, Summary_Statistics.Gap_Statistic(end), 'replicates', 10);
    color = {'k','r','b','m','g','c','y', [0.5, 0, 0.5], [1, 0.5, 0], [0.5, 0.5, 0.5]};
    for jj = 1:Summary_Statistics.Gap_Statistic(end)
        viscircles(pos(idx == jj, 1:2), radii(idx==jj), 'linewidth', 0.75, 'color', color{jj});
    end
   xlim([0 500]);
   ylim([-250, 250]);
   box on; axis square;
   xlabel('$x$ ($\mu$m)', 'interpreter', 'latex', 'fontsize', 16);
   ylabel('$y$ ($\mu$m)', 'interpreter', 'latex', 'fontsize', 16);
   title('Check Gap Statistic', 'interpreter', 'latex', 'fontsize', 16);
   saveas(gcf, [folder_location,'p-Gap_Statistic_Check.svg'], 'svg');
   
   figure('units', 'inches', 'position', [0,0,7,5]);
   hold on;
    idx = kmeans(pos, Summary_Statistics.Gap_Statistic_Standard_1_Error(end));
    for jj = 1:Summary_Statistics.Gap_Statistic_Standard_1_Error(end)
        viscircles(pos(idx == jj, 1:2), radii(idx==jj), 'linewidth', 0.75, 'color', color{jj});
    end
   xlim([0 500]);
   ylim([-250, 250]);
   box on; axis square;
   xlabel('$x$ ($\mu$m)', 'interpreter', 'latex', 'fontsize', 16);
   ylabel('$y$ ($\mu$m)', 'interpreter', 'latex', 'fontsize', 16);
   title('Check Gap Statistic, Standard 1 Error', 'interpreter', 'latex', 'fontsize', 16);
   saveas(gcf, [folder_location,'q-Gap_Statistic_Standard_1_Error_Check.svg'], 'svg');
end

if isfield(Summary_Statistics, 'Porosity_Measure')
    figure('units', 'inches', 'position', [0,0,7,5]);
    p1 = plot((1:length(Summary_Statistics.Porosity_Measure))', Summary_Statistics.Porosity_Measure, 'ob');
   hold on;
   p2 = fplot(@(x) 0*x + mean(Summary_Statistics.Porosity_Measure), [1, length(Summary_Statistics.Porosity_Measure)], '--k');
    legend([p1(1), p2(1)], 'Statistic for Particular Realization', 'Mean for all Realizations');
   xlabel('Realization', 'interpreter', 'latex', 'fontsize', 16);
   ylabel('Measure, ($\mu$m)', 'interpreter', 'latex', 'fontsize', 16);
   title('Average Porosity Measure (Max Y-Width/No. of Cells)', 'interpreter', 'latex', 'fontsize', 16);
   saveas(gcf, [folder_location,'r-Porosity_Measure.svg'], 'svg');
end

if isfield(Summary_Statistics, 'Percent_Target_Area_Covered')
    figure('units', 'inches', 'position', [0,0,7,5]);
    p1 = plot((1:length(Summary_Statistics.Percent_Target_Area_Covered))', 100.*Summary_Statistics.Percent_Target_Area_Covered, 'ob');
   hold on;
   p2 = fplot(@(x) 0*x + 100.*mean(Summary_Statistics.Percent_Target_Area_Covered), [1, length(Summary_Statistics.Percent_Target_Area_Covered)], '--k');
    legend([p1(1), p2(1)], 'Statistic for Particular Realization', 'Mean for all Realizations');
   xlabel('Realization', 'interpreter', 'latex', 'fontsize', 16);
   ylabel('Percent covered', 'interpreter', 'latex', 'fontsize', 16);
   title('Average Percent of Target Corridor Covered', 'interpreter', 'latex', 'fontsize', 16);
   saveas(gcf, [folder_location,'s-Percent_Target_Corridor_Covered.svg'], 'svg');
   
   figure('units', 'inches', 'position', [0,0,7,5]);
   p1 = plot((1:length(Summary_Statistics.Percent_Cell_Area_Outside_Target_Sites))', 100.*Summary_Statistics.Percent_Cell_Area_Outside_Target_Sites, 'ob');
   hold on;
   p2 = fplot(@(x) 0*x + 100.*mean(Summary_Statistics.Percent_Cell_Area_Outside_Target_Sites), [1, length(Summary_Statistics.Percent_Cell_Area_Outside_Target_Sites)], '--k');
   legend([p1(1), p2(1)], 'Statistic for Particular Realization', 'Mean for all Realizations');
   xlabel('Realization', 'interpreter', 'latex', 'fontsize', 16);
   ylabel('Percent', 'interpreter', 'latex', 'fontsize', 16);
   title('Average Percent of Cell Area outside Target Corridor', 'interpreter', 'latex', 'fontsize', 16);
   saveas(gcf, [folder_location,'t-Percent_Cell_Area_Outside_Target_Corridor.svg'], 'svg');
end

if isfield(Summary_Statistics, 'Average_Min_Distance_from_Follower_to_Leader_per_Realization')
    figure('units', 'inches', 'position', [0,0,7,5]);
    p1 = plot((1:length(Summary_Statistics.Average_Min_Distance_from_Follower_to_Leader_per_Realization))', Summary_Statistics.Average_Min_Distance_from_Follower_to_Leader_per_Realization, 'ob');
   hold on;
   p2 = fplot(@(x) 0*x + mean(Summary_Statistics.Average_Min_Distance_from_Follower_to_Leader_per_Realization), [1, length(Summary_Statistics.Average_Min_Distance_from_Follower_to_Leader_per_Realization)], '--k');
    legend([p1(1), p2(1)], 'Statistic for Particular Realization', 'Mean for all Realizations');
   xlabel('Realization', 'interpreter', 'latex', 'fontsize', 16);
   ylabel('Distance ($\mu$m)', 'interpreter', 'latex', 'fontsize', 16);
   title('Average Nearest Dist from Follower to Leader', 'interpreter', 'latex', 'fontsize', 16);
   saveas(gcf, [folder_location,'u-Ave_Nearest_Dist_from_Follower_to_Leader.svg'], 'svg');
   
   figure('units', 'inches', 'position', [0,0,7,5]);
    p1 = plot((1:length(Summary_Statistics.Average_Min_Distance_from_Leader_to_Follower_per_Realization))', Summary_Statistics.Average_Min_Distance_from_Leader_to_Follower_per_Realization, 'ob');
   hold on;
   p2 = fplot(@(x) 0*x + mean(Summary_Statistics.Average_Min_Distance_from_Leader_to_Follower_per_Realization), [1, length(Summary_Statistics.Average_Min_Distance_from_Leader_to_Follower_per_Realization)], '--k');
    legend([p1(1), p2(1)], 'Statistic for Particular Realization', 'Mean for all Realizations');
   xlabel('Realization', 'interpreter', 'latex', 'fontsize', 16);
   ylabel('Distance ($\mu$m)', 'interpreter', 'latex', 'fontsize', 16);
   title('Average Nearest Dist from Leader to Follower', 'interpreter', 'latex', 'fontsize', 16);
   saveas(gcf, [folder_location,'uu-Ave_Nearest_Dist_from_Leader_to_Follower.svg'], 'svg');
   
   figure('units', 'inches', 'position', [0,0,7,5]);
    p1 = plot((1:length(Summary_Statistics.Average_Min_Distance_from_Leader_to_Leader_per_Realization))', Summary_Statistics.Average_Min_Distance_from_Leader_to_Leader_per_Realization, 'ob');
   hold on;
   p2 = fplot(@(x) 0*x + mean(Summary_Statistics.Average_Min_Distance_from_Leader_to_Leader_per_Realization), [1, length(Summary_Statistics.Average_Min_Distance_from_Leader_to_Leader_per_Realization)], '--k');
    legend([p1(1), p2(1)], 'Statistic for Particular Realization', 'Mean for all Realizations');
   xlabel('Realization', 'interpreter', 'latex', 'fontsize', 16);
   ylabel('Distance ($\mu$m)', 'interpreter', 'latex', 'fontsize', 16);
   title('Average Nearest Dist from Leader to Another Leader', 'interpreter', 'latex', 'fontsize', 16);
   saveas(gcf, [folder_location,'v-Ave_Nearest_Dist_from_Leader_to_Leader.svg'], 'svg');
end 

if isfield(Summary_Statistics, 'Number_of_Single_Leader_Streams')
    figure('units', 'inches', 'position', [0,0,7,5]);
    p1 = plot((1:length(Summary_Statistics.Number_of_Single_Leader_Streams))', Summary_Statistics.Number_of_Single_Leader_Streams, 'ob');
   hold on;
   p2 = fplot(@(x) 0*x + mean(Summary_Statistics.Number_of_Single_Leader_Streams), [1, length(Summary_Statistics.Number_of_Single_Leader_Streams)], '--k');
    legend([p1(1), p2(1)], 'Statistic for Particular Realization', 'Mean for all Realizations');
   xlabel('Realization', 'interpreter', 'latex', 'fontsize', 16);
   ylabel('Number of Trails', 'interpreter', 'latex', 'fontsize', 16);
   title('Number of Trails with a Single Leader', 'interpreter', 'latex', 'fontsize', 16);
   saveas(gcf, [folder_location,'w-Number_of_Single_Leader_Streams.svg'], 'svg');
   
   figure('units', 'inches', 'position', [0,0,7,5]);
    p1 = plot((1:length(Summary_Statistics.Number_of_Single_Leader_Streams_Filtered_by_Density_Stat))', Summary_Statistics.Number_of_Single_Leader_Streams_Filtered_by_Density_Stat, 'ob');
   hold on;
   p2 = fplot(@(x) 0*x + mean(Summary_Statistics.Number_of_Single_Leader_Streams_Filtered_by_Density_Stat), [1, length(Summary_Statistics.Number_of_Single_Leader_Streams_Filtered_by_Density_Stat)], '--k');
    legend([p1(1), p2(1)], 'Statistic for Particular Realization', 'Mean for all Realizations');
   xlabel('Realization', 'interpreter', 'latex', 'fontsize', 16);
   ylabel('Number of Trails', 'interpreter', 'latex', 'fontsize', 16);
   title('Number of Trails with a One Leader, Filtered by Density', 'interpreter', 'latex', 'fontsize', 16);
   saveas(gcf, [folder_location,'x-Number_of_Single_Leader_Streams_Filtered_by_Density.svg'], 'svg');
end

if isfield(Summary_Statistics, 'Number_of_Outlier_Cells')
   figure('units', 'inches', 'position', [0,0,7,5]);
   p1 = plot((1:length(Summary_Statistics.Number_of_Outlier_Cells))', Summary_Statistics.Number_of_Outlier_Cells, 'ob');
   hold on;
   p2 = fplot(@(x) 0*x + mean(Summary_Statistics.Number_of_Outlier_Cells), [1, length(Summary_Statistics.Number_of_Outlier_Cells)], '--k');
   legend([p1(1), p2(1)], 'Statistic for Particular Realization', 'Mean for all Realizations');
   xlabel('Realization', 'interpreter', 'latex', 'fontsize', 16);
   ylabel('Number of Trails', 'interpreter', 'latex', 'fontsize', 16);
   title('Number of Cells outside the Target Corridor', 'interpreter', 'latex', 'fontsize', 16);
   saveas(gcf, [folder_location,'y-Number_of_Outlier_Cells.svg'], 'svg');
end

if isfield(Summary_Statistics, 'Max_Leader_X_Position')
   figure('units', 'inches', 'position', [0,0,7,5]);
   p1 = plot((1:length(Summary_Statistics.Max_Leader_X_Position))', Summary_Statistics.Max_Leader_X_Position, 'ob');
   hold on;
   p2 = fplot(@(x) 0*x + mean(Summary_Statistics.Max_Leader_X_Position), [1, length(Summary_Statistics.Max_Leader_X_Position)], '--k');
   legend([p1(1), p2(1)], 'Statistic for Particular Realization', 'Mean for all Realizations');
   xlabel('Realization', 'interpreter', 'latex', 'fontsize', 16);
   ylabel('Distance ($\mu$m)', 'interpreter', 'latex', 'fontsize', 16);
   title('Max Leader $x$-position', 'interpreter', 'latex', 'fontsize', 16);
   saveas(gcf, [folder_location,'z-Max_Leader_X_Position.svg'], 'svg');
   
   figure('units', 'inches', 'position', [0,0,7,5]);
   p1 = plot((1:length(Summary_Statistics.Max_Follower_X_Position))', Summary_Statistics.Max_Follower_X_Position, 'ob');
   hold on;
   p2 = fplot(@(x) 0*x + mean(Summary_Statistics.Max_Follower_X_Position), [1, length(Summary_Statistics.Max_Follower_X_Position)], '--k');
   legend([p1(1), p2(1)], 'Statistic for Particular Realization', 'Mean for all Realizations');
   xlabel('Realization', 'interpreter', 'latex', 'fontsize', 16);
   ylabel('Distance ($\mu$m)', 'interpreter', 'latex', 'fontsize', 16);
   title('Max Follower $x$-position', 'interpreter', 'latex', 'fontsize', 16);
   saveas(gcf, [folder_location,'aa-Max_Follower_X_Position.svg'], 'svg');
   
   figure('units', 'inches', 'position', [0,0,7,5]);
   p1 = plot((1:length(Summary_Statistics.Average_Leader_Speed_per_Realization))', Summary_Statistics.Average_Leader_Speed_per_Realization, 'ob');
   hold on;
   p2 = fplot(@(x) 0*x + mean(Summary_Statistics.Average_Leader_Speed_per_Realization), [1, length(Summary_Statistics.Average_Leader_Speed_per_Realization)], '--k');
   legend([p1(1), p2(1)], 'Statistic for Particular Realization', 'Mean for all Realizations');
   xlabel('Realization', 'interpreter', 'latex', 'fontsize', 16);
   ylabel('Speed ($\mu$m/min)', 'interpreter', 'latex', 'fontsize', 16);
   title('Average Leader Speed', 'interpreter', 'latex', 'fontsize', 16);
   saveas(gcf, [folder_location,'ab-Average_Leader_Speed_per_Realization.svg'], 'svg');
   
   figure('units', 'inches', 'position', [0,0,7,5]);
   p1 = plot((1:length(Summary_Statistics.Average_Follower_Speed_per_Realization))', Summary_Statistics.Average_Follower_Speed_per_Realization, 'ob');
   hold on;
   p2 = fplot(@(x) 0*x + mean(Summary_Statistics.Average_Follower_Speed_per_Realization), [1, length(Summary_Statistics.Average_Follower_Speed_per_Realization)], '--k');
   legend([p1(1), p2(1)], 'Statistic for Particular Realization', 'Mean for all Realizations');
   xlabel('Realization', 'interpreter', 'latex', 'fontsize', 16);
   ylabel('Speed ($\mu$m/min)', 'interpreter', 'latex', 'fontsize', 16);
   title('Average Follower Speed', 'interpreter', 'latex', 'fontsize', 16);
   saveas(gcf, [folder_location,'ac-Average_Follower_Speed_per_Realization.svg'], 'svg');
   
   figure('units', 'inches', 'position', [0,0,7,5]);
   p1 = plot((1:length(Summary_Statistics.Min_Leader_X_Position))', Summary_Statistics.Min_Leader_X_Position, 'ob');
   hold on;
   p2 = fplot(@(x) 0*x + mean(Summary_Statistics.Min_Leader_X_Position), [1, length(Summary_Statistics.Min_Leader_X_Position)], '--k');
   legend([p1(1), p2(1)], 'Statistic for Particular Realization', 'Mean for all Realizations');
   xlabel('Realization', 'interpreter', 'latex', 'fontsize', 16);
   ylabel('Distance ($\mu$m)', 'interpreter', 'latex', 'fontsize', 16);
   title('Min Leader $x$-Position', 'interpreter', 'latex', 'fontsize', 16);
   saveas(gcf, [folder_location,'ad-Min_Leader_X_Position.svg'], 'svg');
end