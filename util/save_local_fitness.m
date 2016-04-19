%% plot local fitness and save them to disk

function [] = save_local_fitness(ct_range)

global output_only_cell_types occupation mid_z fitness_1 cell_type_color_map tick num_iters path;

% cell type local fitness (2-D slice)
for ct = ct_range
    
    % process only designated cell types
    if ~numel(find(output_only_cell_types == ct))
        continue;
    end
    
    % set up an invisible figure
    h = figure('Visible', 'off');
    
    % plot the image
    imagesc((occupation(:,:,mid_z) == ct) .* fitness_1(:,:,mid_z));
    colormap(grade(cell_type_color_map(ct,:), 256));
    cbfreeze(colorbar);
    freezeColors;
    set(gca, 'YDir', 'normal');
    axis square;
    title({['Local Fitness of Cell Type ' num2str(ct)] ['t=' num2str(tick)]});
    
    % create figure file name and save figure
    len_num_iters = length(num2str(num_iters));
    len_tick = length(num2str(tick));
    len_pad = len_num_iters - len_tick;
    pad = '';
    for pi = 1 : len_pad
        pad = strcat('0', pad);
    end
    fig_file = strcat(path, '/', 'local_fitness_', num2str(ct), '_', pad, num2str(tick));
    saveas(gcf, fig_file, 'pdf');
    
    % delete figure
    close(h);
    
end

end  % function save_local_fitness

