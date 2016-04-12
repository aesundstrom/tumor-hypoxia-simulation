%% plot neighborhood fitness and save them to disk
function [] = save_neighborhood_fitness(ct_range)

global output_only_cell_types occupation fitness_n cell_type_color_map tick num_iters path;

% cell type neighborhood fitness (2-D slice)
for ct = ct_range
    
    % process only designated cell types
    if ~numel(find(output_only_cell_types == ct))
        continue;
    end
    
    % set up an invisible figure
    h = figure('Visible', 'off');
    
    % plot the image
    imagesc((occupation(:,:,mid_z) > 1) .* fitness_n(:,:,mid_z,ct));
    colormap(grade(cell_type_color_map(ct,:), 256));
    cbfreeze(colorbar);
    freezeColors;
    set(gca, 'YDir', 'normal');
    axis square;
    title({['Neighborhood Fitness of Cell Type ' num2str(ct)] ['t=' num2str(tick)]});
    
    % create figure file name and save figure
    len_num_iters = length(num2str(num_iters));
    len_tick = length(num2str(tick));
    len_pad = len_num_iters - len_tick;
    pad = '';
    for pi = 1 : len_pad
        pad = strcat('0', pad);
    end
    fig_file = strcat(path, '/', 'neighborhood_fitness_', num2str(ct), '_', pad, num2str(tick));
    saveas(gcf, fig_file, 'pdf');
    
    % delete figure
    close(h);
    
end

end  % function save_neighborhood_fitness

