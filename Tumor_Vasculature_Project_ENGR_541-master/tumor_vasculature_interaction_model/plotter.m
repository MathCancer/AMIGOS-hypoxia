hold off
clear all
%load vasculature, plot with contour
%AXIS IMAGE
for i=1:40
    filename = sprintf('output%08u_vasculature.mat', i);
    file  = sprintf('output%08u.xml', i);
    VASC = read_microenvironment(filename);
    MCDS = read_MultiCellDS_xml(file); % 08u, i 
    O2 =  MCDS.continuum_variables(1).data;
    VEGF = MCDS.continuum_variables(2).data;
    xcords = MCDS.mesh.X_coordinates;
    ycords = MCDS.mesh.Y_coordinates;
    PLOT = figure;
    set(PLOT, 'Visible', 'off');
%     contourf(VASC.X, VASC.Y, VASC.data{1,1});
    time = string(sprintf('%d', MCDS.metadata.current_time/60));
%     title({'Vasculature' ; 'Time = ' + time + ' hrs'});
%     colorbar;
%     set(colorbar, 'Ylim', [0,1])
%     axis image;
%     hold on;
%     plot(MCDS.discrete_cells.state.position(:,1), ...
%         MCDS.discrete_cells.state.position(:,2), 'ko');
%     vasculature = sprintf('VASC%08u.png', i);
%     saveas(PLOT, vasculature);
%     clf
%     contourf(xcords, ycords, O2)
%     title({'Oxygen' ; 'Time = ' + time + ' hrs'});
%     colorbar;
%     set(colorbar, 'Ylim', [0,38])
%     axis image;
%     hold on;
%     plot(MCDS.discrete_cells.state.position(:,1), ...
%         MCDS.discrete_cells.state.position(:,2), 'ko');
%     O2_file = sprintf('O2_%08u.png', i);
%     saveas(PLOT, O2_file)
    clf
    contourf(xcords, ycords, VEGF)
    title({'VEGF' ; 'Time = ' + time + ' hrs'});
    colorbar;
    set(colorbar, 'Ylim', [0,1])
    axis image;
    hold on;
    plot(MCDS.discrete_cells.state.position(:,1), ...
        MCDS.discrete_cells.state.position(:,2), 'ko');
    VEGF_file = sprintf('VEGF_%08u.png', i);
    saveas(PLOT, VEGF_file)
    clf
end