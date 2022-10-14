function MCDS = plot_cell_ecm_interactions
% close all;
% Find out how many XML files there are
listing = dir('/Users/duncanmartinson/Documents/PhysiCell/output/*.xml');
N = length(listing);
% Load in the MCDS files. The first two files are 'final.xml' and
% 'initial.xml', so start from i = 3:
% Load the first MCDS File
MCDS = read_MultiCellDS_xml(listing(end).name, '/Users/duncanmartinson/Documents/PhysiCell/output');
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

a = VideoWriter('ECM-field-progression', 'MPEG-4');
a.Quality = 100;
a.FrameRate = 5;
open(a);
for i = 3:N
    % Load the MCDS file:
    MCDS = read_MultiCellDS_xml(listing(i).name, '/Users/duncanmartinson/Documents/PhysiCell/output');
    % Extract the cells
    cells = MCDS.discrete_cells;
    Leaders = (cells.metadata.type == 0);
    Followers = (cells.metadata.type == 1);
    ECM = (cells.metadata.type == 2);
    % Plot the ECM vector field and then the cells themselves as a scatter
    % plot:
    quiver(cells.state.position(ECM, 1), cells.state.position(ECM, 2), ...
    cells.phenotype.motility.motility_vector(ECM, 1), cells.phenotype.motility.motility_vector(ECM,2),...
    'k', 'AutoScaleFactor', 0.2);
    hold on;
    plot(cells.state.position(Leaders, 1), cells.state.position(Leaders, 2), 'ok');
    plot(cells.state.position(Followers, 1), cells.state.position(Followers, 2), 'or');
    set(get(ax1, 'title'), 'string', ['$t$ = ',num2str(MCDS.metadata.current_time), ' ', MCDS.metadata.time_units]);
    ax1.Title.Interpreter = 'latex';
    ax1.FontSize = 20;
    xlabel('$x$', 'interpreter', 'latex', 'fontsize', 20);
    ylabel('$y$', 'interpreter', 'latex', 'fontsize', 20);  
    hold off;
    
    xlim([Xmin, Xmax]);
    ylim([Ymin, Ymax]);
    
    frame = getframe(fig1);
    writeVideo(a, frame);
end % for i
close(a);

fig2 = figure;
% Plot ECM vector field
quiver(cells.state.position(ECM, 1), cells.state.position(ECM, 2), ...
    cells.phenotype.motility.motility_vector(ECM, 1), cells.phenotype.motility.motility_vector(ECM,2),...
    'k', 'AutoScaleFactor', 0.25);
hold on;
plot(cells.state.position(Leaders, 1), cells.state.position(Leaders, 2), 'ok');
plot(cells.state.position(Followers, 1), cells.state.position(Followers, 2), 'or');
xlim([Xmin, Xmax]);
ylim([Ymin, Ymax]);
xlabel('$x$', 'interpreter', 'latex', 'fontsize', 20);
ylabel('$y$', 'interpreter', 'latex', 'fontsize', 20);   
title(['$t$ = ',num2str(MCDS.metadata.current_time), ' ', MCDS.metadata.time_units],'interpreter', 'latex', 'fontsize', 20);
end