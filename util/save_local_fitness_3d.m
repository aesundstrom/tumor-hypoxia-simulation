%% plot local fitness 3D and save them to disk
function [] = save_local_fitness_3d(ct_range)

global output_only_cell_types occupation fitness_1 cell_type_color_map x_dim y_dim z_dim tick num_iters path;

% cell type 3D plot
for ct = ct_range
    
    % process only designated cell types
    if ~numel(find(output_only_cell_types == ct))
        continue;
    end
    
    % set up an invisible figure
    h = figure('Visible', 'off');
    
    % plot the image
    occ_ct = occupation == ct;
    fit_ct = fitness_1(occ_ct);
    min_fit_ct = min(fit_ct);
    max_fit_ct = max(fit_ct);
    norm_fit_ct = fit_ct;
    if min_fit_ct < max_fit_ct
        norm_fit_ct = (fit_ct - min_fit_ct) ./ (max_fit_ct - min_fit_ct);
    end
    color_density = 256;
    color_grade = grade(cell_type_color_map(ct,:), color_density);
    color_idx = floor((color_density - 1) * norm_fit_ct) + 1;
    scatter_colors = color_grade(color_idx,:);
    [x,y,z] = ind2sub(size(occ_ct), find(occ_ct));
    % note the axis order for plotting: <y,x,z>
    scatter3(y, x, z, 20, scatter_colors, 'filled', 's');
    colormap(grade(cell_type_color_map(ct,:), color_density));
    drawnow;
    ch = colorbar;
    if min_fit_ct < max_fit_ct
        caxis([min_fit_ct max_fit_ct]);
    end
    cbfreeze(ch);
    xlabel('x');
    ylabel('y');
    zlabel('z');
    if z_dim > 1
        axis([1 x_dim 1 y_dim 1 z_dim]);
    else
        axis([1 x_dim 1 y_dim]);
    end
    axis square;
    title({['Local Fitness of Cell Type ' num2str(ct) ' (3D View)'] ['t=' num2str(tick)]});
    
    % create figure file name and save figure
    len_num_iters = length(num2str(num_iters));
    len_tick = length(num2str(tick));
    len_pad = len_num_iters - len_tick;
    pad = '';
    for pi = 1 : len_pad
        pad = strcat('0', pad);
    end
    fig_file = strcat(path, '/', 'local_fitness_3d_', num2str(ct), '_', pad, num2str(tick));
    saveas(gcf, fig_file, 'pdf');
    
    % delete figure
    close(h);
    
end

end  % function save_local_fitness_3d

