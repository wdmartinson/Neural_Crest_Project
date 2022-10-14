function MCDS = write_vicsek_videos
% Find out how many XML files there are
listing = dir('/Users/duncanmartinson/Documents/PhysiCell/output/*.xml');
N = length(listing);
% Load in the MCDS files. The first two outputs are 'final.xml' and
% 'initial.xml', so start from i = 3:
% Load the first MCDS File
MCDS = read_MultiCellDS_xml(listing(2).name, '/Users/duncanmartinson/Documents/PhysiCell/output');
dx = MCDS.mesh.X_coordinates(2)-MCDS.mesh.X_coordinates(1);
dy = MCDS.mesh.Y_coordinates(2)-MCDS.mesh.Y_coordinates(1);
Xmin = MCDS.mesh.X_coordinates(1) - dx/2;
Ymin = MCDS.mesh.Y_coordinates(1) - dy/2;
Xmax = MCDS.mesh.X_coordinates(end) + dx/2;
Ymax = MCDS.mesh.Y_coordinates(end) + dy/2;

fig1 = figure;
    graph1 = plot(NaN);
    xlim([Xmin, Xmax]);
    ylim([Ymin, Ymax]);
    xlabel('$x$', 'interpreter', 'latex', 'fontsize', 20);
    ylabel('$y$', 'interpreter', 'latex', 'fontsize', 20);   
    axis tight;
    set(gca, 'nextplot', 'replacechildren', 'visible','on');
    ax1 = gca;
    title('','interpreter', 'latex', 'fontsize', 20);
    set(get(ax1, 'title'), 'string', ['$t$ = ',num2str(MCDS.metadata.current_time), ' ', MCDS.metadata.time_units]);

a = VideoWriter('vicsek_model_test', 'MPEG-4');
a.Quality = 100;
a.FrameRate = 5;
open(a);
for i = 3:N
    % Load the MCDS file:
    MCDS = read_MultiCellDS_xml(listing(i).name, '/Users/duncanmartinson/Documents/PhysiCell/output');
    % Extract the cells
    cells = MCDS.discrete_cells;
    % Plot the cells as a scatter plot
    plot(ax1, cells.state.position(:,1), cells.state.position(:,2), 'ok', 'MarkerFaceColor','black');
    hold on;
    quiver(cells.state.position(:, 1), cells.state.position(:, 2), ...
        cells.phenotype.motility.motility_vector(:, 1), cells.phenotype.motility.motility_vector(:,2),...
        'k', 'AutoScaleFactor', 0.2);
    set(get(ax1, 'title'), 'string', ['$t$ = ',num2str(MCDS.metadata.current_time), ' ', MCDS.metadata.time_units]);
    hold off;
    
    xlim([Xmin, Xmax]);
    ylim([Ymin, Ymax]);
    
    frame = getframe(fig1);
    writeVideo(a, frame);
end % for i
close(a);

fig2 = figure;
plot(cells.state.position(:,1), cells.state.position(:,2), 'ok', 'MarkerFaceColor','black');
xlim([Xmin, Xmax]);
ylim([Ymin, Ymax]);
xlabel('$x$', 'interpreter', 'latex', 'fontsize', 20);
ylabel('$y$', 'interpreter', 'latex', 'fontsize', 20);   
axis tight;
title(['$t$ = ',num2str(MCDS.metadata.current_time), ' ', MCDS.metadata.time_units],'interpreter', 'latex', 'fontsize', 20);
hold on;
quiver(MCDS.discrete_cells.state.position(:, 1), MCDS.discrete_cells.state.position(:, 2), ...
    MCDS.discrete_cells.phenotype.motility.motility_vector(:, 1), MCDS.discrete_cells.phenotype.motility.motility_vector(:,2),...
    'k', 'AutoScaleFactor', 0.2);
xlim([Xmin, Xmax]);
ylim([Ymin, Ymax]);
hold off;
end