function quick_plot_NCC_experiment
folder_location = './output/';
% folder_location = '../temp2/';

% Find out how many XML files there are
listing = dir(strcat(folder_location,'*.xml'));
N = length(listing);

% Load in the MCDS files. The first three files are 'config.xml', 'final.xml', and
% 'initial.xml', so start from i = 4:
% Load the first MCDS File
MCDS = read_MultiCellDS_xml(listing(end).name, folder_location);
dx = MCDS.mesh.X_coordinates(2)-MCDS.mesh.X_coordinates(1);
dy = MCDS.mesh.Y_coordinates(2)-MCDS.mesh.Y_coordinates(1);
Xmin = MCDS.mesh.X_coordinates(1) - dx/2;
Ymin = MCDS.mesh.Y_coordinates(1) - dy/2;
Xmax = MCDS.mesh.X_coordinates(end) + dx/2;
Ymax = MCDS.mesh.Y_coordinates(end) + dy/2;

% MCDS = read_MultiCellDS_xml(listing(3).name, folder_location);
MCDS = read_MultiCellDS_xml(listing(end).name, folder_location);
% MCDS = read_MultiCellDS_xml(listing(4 + 36*4).name, folder_location);
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
%     caxis([0 FN_max]);
%     set(get(cbar, 'title'), 'string', 'Fibronectin, $\mu$M', 'interpreter', 'latex', 'fontsize', 16);
    ax1.Title.Interpreter = 'latex';
    ax1.FontSize = 16;
    xlim([Xmin, Xmax]);
    ylim([Ymin, Ymax]);
    
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
%     cbar = colorbar; caxis([0 FN_max]);
    hold on;
%     viscircles([cells.state.position(Leaders,1) cells.state.position(Leaders,2)], leader_radii, 'color', 'k', 'linewidth', 0.5);
%     viscircles([cells.state.position(Followers,1) cells.state.position(Followers,2)], follower_radii, 'color', 'r', 'linewidth', 0.5);
% 
%     % Also draw unit vectors to indicate current velocity of the cell:
%     quiver(cells.state.position(Leaders, 1), cells.state.position(Leaders, 2), ...
%         cells.custom.total_velocity(Leaders, 1)./sqrt(sum(cells.custom.total_velocity(Leaders, :).^2, 2)), cells.custom.total_velocity(Leaders,2)./sqrt(sum(cells.custom.total_velocity(Leaders, :).^2, 2)),...
%         'k', 'AutoScaleFactor', 0.2);
%     quiver(cells.state.position(Followers, 1), cells.state.position(Followers, 2), ...
%         cells.custom.total_velocity(Followers, 1)./sqrt(sum(cells.custom.total_velocity(Followers, :).^2, 2)), cells.custom.total_velocity(Followers,2)./sqrt(sum(cells.custom.total_velocity(Followers, :).^2, 2)),...
%         'r', 'AutoScaleFactor', 0.2);
    h1 = viscircles([cells.state.position(Leaders,1) cells.state.position(Leaders,2)], leader_radii, 'color', 'k', 'linewidth', 0.75);
    x_leaders = h1.Children(1).XData;
    y_leaders = h1.Children(1).YData;
    idx = [0, find(isnan(x_leaders))];
    for ii = 1:length(idx)-1
       fill(x_leaders(idx(ii)+1:idx(ii+1)-1), y_leaders(idx(ii)+1:idx(ii+1)-1), 'k', 'facealpha', 0.5);
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
        0, 'color','k');
    quiver(cells.state.position(Followers, 1), cells.state.position(Followers, 2), ...
        15*cells.custom.total_velocity(Followers, 1)./sqrt(sum(cells.custom.total_velocity(Followers, :).^2, 2)), 15*cells.custom.total_velocity(Followers,2)./sqrt(sum(cells.custom.total_velocity(Followers, :).^2, 2)),...
        0, 'color','r');
     quiver(cells.state.position(Puncta1, 1), cells.state.position(Puncta1, 2), ...
        15*cos(cells.custom.time_until_next_filopodia_drop(Puncta1)'), 15*sin(cells.custom.time_until_next_filopodia_drop(Puncta1)'),...
        0, 'color', 'b');
    % Housekeeping: Keeping the video quality the same throughout
    xlim([Xmin, Xmax]); ylim([Ymin, Ymax]);
    set(get(ax1, 'title'), 'string', ['$t$ = ',num2str(MCDS.metadata.current_time), ' ', MCDS.metadata.time_units]);
%     set(get(cbar, 'title'), 'string', 'Fibronectin, $\mu$M', 'interpreter', 'latex', 'fontsize', 16);
    ax1.Title.Interpreter = 'latex';
    ax1.FontSize = 20;
    set(ax1, 'TickLabelInterpreter','latex','fontsize',20);
    xlabel('$x$, $\mu$m', 'interpreter', 'latex', 'fontsize', 20);
    ylabel('$y$, $\mu$m', 'interpreter', 'latex', 'fontsize', 20); 
end