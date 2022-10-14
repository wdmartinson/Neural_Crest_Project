function plot_ncc_model
% close all;
% Find out how many XML files there are
listing = dir('/Users/duncanmartinson/Documents/PhysiCell/output/*.xml');
N = length(listing);

Video_title = '106-leader_DAN-VEGF_migration_rule_with_cell_cell_repulsion_DAN_is_present_DAN_width_150_DAN_uptake_25e-3_VEGF_uptake_1e-3_new_cells_can_enter_cell_entrance_width_200_cell_IC_along_entrance_width';
% Load in the MCDS files. The first two files are 'final.xml' and
% 'initial.xml', so start from i = 3:
% Load the first MCDS File
MCDS = read_MultiCellDS_xml(listing(2).name, '/Users/duncanmartinson/Documents/PhysiCell/output');
dx = MCDS.mesh.X_coordinates(2)-MCDS.mesh.X_coordinates(1);
dy = MCDS.mesh.Y_coordinates(2)-MCDS.mesh.Y_coordinates(1);
Xmin = MCDS.mesh.X_coordinates(1) - dx/2;
Ymin = MCDS.mesh.Y_coordinates(1) - dy/2;
Xmax = MCDS.mesh.X_coordinates(end) + dx/2;
Ymax = MCDS.mesh.Y_coordinates(end) + dy/2;
VEGF_max_init = max(max(MCDS.continuum_variables(2).data));
if length(MCDS.continuum_variables)>2
    DAN_max_init = max(max(MCDS.continuum_variables(3).data));
end % if there's DAN
MCDS = read_MultiCellDS_xml(listing(end).name, '/Users/duncanmartinson/Documents/PhysiCell/output');
FN_max = max(max(MCDS.continuum_variables(1).data));

fig1 = figure;
    imagesc('xdata', MCDS.mesh.X_coordinates, 'ydata', fliplr(MCDS.mesh.Y_coordinates), 'cdata', NaN);
    xlim([Xmin, Xmax]);
    ylim([Ymin, Ymax]);
    xlabel('$x$, $\mu$m', 'interpreter', 'latex', 'fontsize', 20);
    ylabel('$y$, $\mu$m', 'interpreter', 'latex', 'fontsize', 20);
    
    axis tight;
    set(gca, 'nextplot', 'replacechildren', 'visible','on');
    ax1 = gca;
    title('','interpreter', 'latex', 'fontsize', 20);
    set(get(ax1, 'title'), 'string', ['$t$ = ',num2str(MCDS.metadata.current_time), ' ', MCDS.metadata.time_units]);
    cbar = colorbar;
    caxis([0 FN_max]);
    set(get(cbar, 'title'), 'string', 'Fibronectin, $\mu$M', 'interpreter', 'latex', 'fontsize', 16);
%     cla;
figure_properties = get(gcf);
% Calculate total marker size so that it has same aspect ratio as domain:
% Total area of pixels:
pixel_area = prod(figure_properties.InnerPosition(3:4));
domain_area = (Xmax-Xmin)*(Ymax-Ymin);
% Domain area/pixel_area = cell_radius/Marker_area
a = VideoWriter([Video_title,'-FN_field'], 'MPEG-4');
a.Quality = 100;
a.FrameRate = 5;
open(a);
for i = 3:N
    cla;
    % Load the MCDS file:
    MCDS = read_MultiCellDS_xml(listing(i).name, '/Users/duncanmartinson/Documents/PhysiCell/output');
    % Extract the cells
    cells = MCDS.discrete_cells;
    Leaders = (cells.metadata.type == 0);
    Followers = (cells.metadata.type == 1);
    % Extract the microenvironment:
    FN_field = MCDS.continuum_variables(1).data;
    % Plot the FN field and then the cells themselves as a scatter
    % plot:
    Leader_marker_sizes = (domain_area/pixel_area)^(-1)*2*cells.phenotype.geometrical_properties.lengths.radius(Leaders).^2;
    Follower_marker_sizes = (domain_area/pixel_area)^(-1)*2*cells.phenotype.geometrical_properties.lengths.radius(Followers).^2;
%     Leader_marker_sizes = (domain_area/pixel_area)^(-1)*100*ones(size(cells.phenotype.geometrical_properties.lengths.radius(Leaders)));
%     Follower_marker_sizes = (domain_area/pixel_area)^(-1)*100*ones(size(cells.phenotype.geometrical_properties.lengths.radius(Followers)));
    imagesc('xdata', MCDS.mesh.X_coordinates, 'ydata', fliplr(MCDS.mesh.Y_coordinates), 'cdata', flipud(FN_field));
    cbar = colorbar; caxis([0 FN_max]);
    hold on;
    scatter(cells.state.position(Leaders, 1), cells.state.position(Leaders, 2), Leader_marker_sizes, 'ok');
    scatter(cells.state.position(Followers, 1), cells.state.position(Followers, 2), Follower_marker_sizes, 'or');
    % Also draw unit vectors to indicate current velocity of the cell:
    quiver(cells.state.position(Leaders, 1), cells.state.position(Leaders, 2), ...
        cells.custom.total_velocity(Leaders, 1)./sqrt(sum(cells.custom.total_velocity(Leaders, :).^2, 2)), cells.custom.total_velocity(Leaders,2)./sqrt(sum(cells.custom.total_velocity(Leaders, :).^2, 2)),...
        'k', 'AutoScaleFactor', 0.2);
    quiver(cells.state.position(Followers, 1), cells.state.position(Followers, 2), ...
        cells.custom.total_velocity(Followers, 1)./sqrt(sum(cells.custom.total_velocity(Followers, :).^2, 2)), cells.custom.total_velocity(Followers,2)./sqrt(sum(cells.custom.total_velocity(Followers, :).^2, 2)),...
        'r', 'AutoScaleFactor', 0.2);
    % Housekeeping: Keeping the video quality the same throughout
    xlim([Xmin, Xmax]); ylim([Ymin, Ymax]);
    set(get(ax1, 'title'), 'string', ['$t$ = ',num2str(MCDS.metadata.current_time), ' ', MCDS.metadata.time_units]);
    set(get(cbar, 'title'), 'string', 'Fibronectin, $\mu$M', 'interpreter', 'latex', 'fontsize', 16);
    ax1.Title.Interpreter = 'latex';
%     ax1.FontSize = 20;
    xlabel('$x$, $\mu$m', 'interpreter', 'latex', 'fontsize', 20);
    ylabel('$y$, $\mu$m', 'interpreter', 'latex', 'fontsize', 20); 
    hold off;
    
    % Take the snapshot
    frame = getframe(fig1);
    writeVideo(a, frame);
end % for i
close(a);

fig2 = figure;
    imagesc('xdata', MCDS.mesh.X_coordinates, 'ydata', fliplr(MCDS.mesh.Y_coordinates), 'cdata', NaN);
    xlim([Xmin, Xmax]);
    ylim([Ymin, Ymax]);
    xlabel('$x$, $\mu$m', 'interpreter', 'latex', 'fontsize', 20);
    ylabel('$y$, $\mu$m', 'interpreter', 'latex', 'fontsize', 20);   
    axis tight;
    set(gca, 'nextplot', 'replacechildren', 'visible','on');
    ax2 = gca;
    title('','interpreter', 'latex', 'fontsize', 20);
    set(get(ax2, 'title'), 'string', ['$t$ = ',num2str(MCDS.metadata.current_time), ' ', MCDS.metadata.time_units]);
    cbar = colorbar;
    caxis([0 VEGF_max_init]);
    set(get(cbar, 'title'), 'string', 'VEGF, $\mu$M', 'interpreter', 'latex', 'fontsize', 16);
 
figure_properties = get(gcf);
% Calculate total marker size so that it has same aspect ratio as domain:
% Total area of pixels:
pixel_area = prod(figure_properties.InnerPosition(3:4));

b = VideoWriter([Video_title,'-VEGF_field'], 'MPEG-4');
b.Quality = 100;
b.FrameRate = 5;
open(b);
for i = 3:N
    cla;
    % Load the MCDS file:
    MCDS = read_MultiCellDS_xml(listing(i).name, '/Users/duncanmartinson/Documents/PhysiCell/output');
    % Extract the cells
    cells = MCDS.discrete_cells;
    Leaders = (cells.metadata.type == 0);
    Followers = (cells.metadata.type == 1);
    % Extract the microenvironment:
    VEGF_field = MCDS.continuum_variables(2).data;
    % Calculate marker sizes for scatter plot:
    Leader_marker_sizes = (domain_area/pixel_area)^(-1)*2*cells.phenotype.geometrical_properties.lengths.radius(Leaders).^2;
    Follower_marker_sizes = (domain_area/pixel_area)^(-1)*2*cells.phenotype.geometrical_properties.lengths.radius(Followers).^2;
%     Leader_marker_sizes = (domain_area/pixel_area)^(-1)*100*ones(size(cells.phenotype.geometrical_properties.lengths.radius(Leaders)));
%     Follower_marker_sizes = (domain_area/pixel_area)^(-1)*100*ones(size(cells.phenotype.geometrical_properties.lengths.radius(Followers)));
    % Plot the FN field and then the cells themselves as a scatter plot:
    imagesc('xdata', MCDS.mesh.X_coordinates, 'ydata', fliplr(MCDS.mesh.Y_coordinates), 'cdata', flipud(VEGF_field));
    cbar = colorbar; caxis([0 VEGF_max_init]);
    hold on;
    scatter(cells.state.position(Leaders, 1), cells.state.position(Leaders, 2), Leader_marker_sizes,'ok');
    scatter(cells.state.position(Followers, 1), cells.state.position(Followers, 2), Follower_marker_sizes, 'or');
    % Also draw unit vectors to indicate current velocity of the cell:
    quiver(cells.state.position(Leaders, 1), cells.state.position(Leaders, 2), ...
        cells.custom.total_velocity(Leaders, 1)./sqrt(sum(cells.custom.total_velocity(Leaders, :).^2, 2)), cells.custom.total_velocity(Leaders,2)./sqrt(sum(cells.custom.total_velocity(Leaders, :).^2, 2)),...
        'k', 'AutoScaleFactor', 0.2);
    quiver(cells.state.position(Followers, 1), cells.state.position(Followers, 2), ...
        cells.custom.total_velocity(Followers, 1)./sqrt(sum(cells.custom.total_velocity(Followers, :).^2, 2)), cells.custom.total_velocity(Followers,2)./sqrt(sum(cells.custom.total_velocity(Followers, :).^2, 2)),...
        'r', 'AutoScaleFactor', 0.2);
    % Housekeeping: Keeping the video quality the same throughout
    xlim([Xmin, Xmax]); ylim([Ymin, Ymax]);
    set(get(ax2, 'title'), 'string', ['$t$ = ',num2str(MCDS.metadata.current_time), ' ', MCDS.metadata.time_units]);
    set(get(cbar, 'title'), 'string', 'VEGF, $\mu$M', 'interpreter', 'latex', 'fontsize', 16);
    ax2.Title.Interpreter = 'latex';
    ax2.FontSize = 20;
    xlabel('$x$, $\mu$m', 'interpreter', 'latex', 'fontsize', 20);
    ylabel('$y$, $\mu$m', 'interpreter', 'latex', 'fontsize', 20); 
    hold off;
    
    % Take the snapshot
    frame = getframe(fig2);
    writeVideo(b, frame);
end % for i
close(b);

% If the MCDS file has a DAN field, retrieve it:
if length(MCDS.continuum_variables)>2
fig3 = figure;
    imagesc('xdata', MCDS.mesh.X_coordinates, 'ydata', fliplr(MCDS.mesh.Y_coordinates), 'cdata', NaN);
    xlim([Xmin, Xmax]);
    ylim([Ymin, Ymax]);
    xlabel('$x$, $\mu$m', 'interpreter', 'latex', 'fontsize', 20);
    ylabel('$y$, $\mu$m', 'interpreter', 'latex', 'fontsize', 20);   
    axis tight;
    set(gca, 'nextplot', 'replacechildren', 'visible','on');
    ax3 = gca;
    title('','interpreter', 'latex', 'fontsize', 20);
    set(get(ax3, 'title'), 'string', ['$t$ = ',num2str(MCDS.metadata.current_time), ' ', MCDS.metadata.time_units]);
    cbar = colorbar;
    caxis([0 DAN_max_init]);
    set(get(cbar, 'title'), 'string', 'DAN, $\mu$M', 'interpreter', 'latex', 'fontsize', 16);
 
figure_properties = get(gcf);
% Calculate total marker size so that it has same aspect ratio as domain:
% Total area of pixels:
pixel_area = prod(figure_properties.InnerPosition(3:4));

c = VideoWriter([Video_title,'-DAN_field'], 'MPEG-4');
c.Quality = 100;
c.FrameRate = 5;
open(c);
for i = 3:N
    cla;
    % Load the MCDS file:
    MCDS = read_MultiCellDS_xml(listing(i).name, '/Users/duncanmartinson/Documents/PhysiCell/output');
    % Extract the cells
    cells = MCDS.discrete_cells;
    Leaders = (cells.metadata.type == 0);
    Followers = (cells.metadata.type == 1);
    % Extract the microenvironment:
    DAN_field = MCDS.continuum_variables(3).data;
    % Calculate marker sizes for scatter plot:
    Leader_marker_sizes = (domain_area/pixel_area)^(-1)*2*cells.phenotype.geometrical_properties.lengths.radius(Leaders).^2;
    Follower_marker_sizes = (domain_area/pixel_area)^(-1)*2*cells.phenotype.geometrical_properties.lengths.radius(Followers).^2;
%     Leader_marker_sizes = (domain_area/pixel_area)^(-1)*100*ones(size(cells.phenotype.geometrical_properties.lengths.radius(Leaders)));
%     Follower_marker_sizes = (domain_area/pixel_area)^(-1)*100*ones(size(cells.phenotype.geometrical_properties.lengths.radius(Followers)));
    % Plot the FN field and then the cells themselves as a scatter plot:
    imagesc('xdata', MCDS.mesh.X_coordinates, 'ydata', fliplr(MCDS.mesh.Y_coordinates), 'cdata', flipud(DAN_field));
    cbar = colorbar; caxis([0 DAN_max_init]);
    hold on;
    scatter(cells.state.position(Leaders, 1), cells.state.position(Leaders, 2), Leader_marker_sizes,'ok');
    scatter(cells.state.position(Followers, 1), cells.state.position(Followers, 2), Follower_marker_sizes, 'or');
    % Also draw unit vectors to indicate current velocity of the cell:
    quiver(cells.state.position(Leaders, 1), cells.state.position(Leaders, 2), ...
        cells.custom.total_velocity(Leaders, 1)./sqrt(sum(cells.custom.total_velocity(Leaders, :).^2, 2)), cells.custom.total_velocity(Leaders,2)./sqrt(sum(cells.custom.total_velocity(Leaders, :).^2, 2)),...
        'k', 'AutoScaleFactor', 0.2);
    quiver(cells.state.position(Followers, 1), cells.state.position(Followers, 2), ...
        cells.custom.total_velocity(Followers, 1)./sqrt(sum(cells.custom.total_velocity(Followers, :).^2, 2)), cells.custom.total_velocity(Followers,2)./sqrt(sum(cells.custom.total_velocity(Followers, :).^2, 2)),...
        'r', 'AutoScaleFactor', 0.2);
    % Housekeeping: Keeping the video quality the same throughout
    xlim([Xmin, Xmax]); ylim([Ymin, Ymax]);
    set(get(ax3, 'title'), 'string', ['$t$ = ',num2str(MCDS.metadata.current_time), ' ', MCDS.metadata.time_units]);
    set(get(cbar, 'title'), 'string', 'DAN, $\mu$M', 'interpreter', 'latex', 'fontsize', 16);
    ax3.Title.Interpreter = 'latex';
%     ax3.FontSize = 20;
    xlabel('$x$, $\mu$m', 'interpreter', 'latex', 'fontsize', 20);
    ylabel('$y$, $\mu$m', 'interpreter', 'latex', 'fontsize', 20); 
    hold off;
    
    % Take the snapshot
    frame = getframe(fig3);
    writeVideo(c, frame);
end % for i
close(c);

end % if there is DAN.

% Simple scatter plot at final time point, with unit vector velocities:
figure;
figure_properties = get(gcf);
% Calculate total marker size so that it has same aspect ratio as domain:
% Total area of pixels:
pixel_area = prod(figure_properties.InnerPosition(3:4));
Leader_marker_sizes = (domain_area/pixel_area)^(-1)*2*cells.phenotype.geometrical_properties.lengths.radius(Leaders).^2;
Follower_marker_sizes = (domain_area/pixel_area)^(-1)*2*cells.phenotype.geometrical_properties.lengths.radius(Followers).^2;
% Leader_marker_sizes = (domain_area/pixel_area)^(-1)*100*ones(size(cells.phenotype.geometrical_properties.lengths.radius(Leaders)));
% Follower_marker_sizes = (domain_area/pixel_area)^(-1)*100*ones(size(cells.phenotype.geometrical_properties.lengths.radius(Followers)));
hold on; box on;
scatter(cells.state.position(Leaders, 1), cells.state.position(Leaders, 2), Leader_marker_sizes,'ok');
scatter(cells.state.position(Followers, 1), cells.state.position(Followers, 2), Follower_marker_sizes, 'or');
quiver(cells.state.position(Leaders, 1), cells.state.position(Leaders, 2), ...
    cells.custom.total_velocity(Leaders, 1)./sqrt(sum(cells.custom.total_velocity(Leaders, :).^2, 2)), cells.custom.total_velocity(Leaders,2)./sqrt(sum(cells.custom.total_velocity(Leaders, :).^2, 2)),...
    'k', 'AutoScaleFactor', 0.2);
quiver(cells.state.position(Followers, 1), cells.state.position(Followers, 2), ...
    cells.custom.total_velocity(Followers, 1)./sqrt(sum(cells.custom.total_velocity(Followers, :).^2, 2)), cells.custom.total_velocity(Followers,2)./sqrt(sum(cells.custom.total_velocity(Followers, :).^2, 2)),...
    'r', 'AutoScaleFactor', 0.2);
xlim([Xmin, Xmax]);
ylim([Ymin, Ymax]);
xlabel('$x$, $\mu$m', 'interpreter', 'latex', 'fontsize', 20);
ylabel('$y$, $\mu$m', 'interpreter', 'latex', 'fontsize', 20);   
title(['$t$ = ',num2str(MCDS.metadata.current_time), ' ', MCDS.metadata.time_units],'interpreter', 'latex', 'fontsize', 20);
end