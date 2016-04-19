%% plot particle concentration and save them to disk

function [] = save_particle_concentration(gray_256, pt_range)

global x_dim y_dim path tick;

% particle type concentrations (2-D slice)
for pt = pt_range
    
    % process only designated particle types
    if ~numel(find(output_only_particle_types == ct))
        continue;
    end
    
    % set up an invisible figure
    h = figure('Visible', 'off');
    
    % plot the image
    con_pt = concentration(:,:,mid_z,pt);
    min_con_pt = min(min(con_pt));
    max_con_pt = max(max(con_pt));
    mesh(con_pt);
    colormap(gray_256);
    min_rel = '=';
    if min_con_pt > 1000
        min_rel = '>';
        min_con_pt = 1000;
    end
    max_rel = '=';
    if max_con_pt > 1000
        max_rel = '>';
        max_con_pt = 1000;
    end
    axis square;
    xlim([1 x_dim]);
    ylim([1 y_dim]);
    title({['Concentration of Particle Type ' num2str(pt)] ['t=' num2str(tick)] ['min' min_rel sprintf('%3.3f', min_con_pt) ', max' max_rel sprintf('%3.3f', max_con_pt)]});
    
    % create figure file name and save figure
    len_num_iters = length(num2str(num_iters));
    len_tick = length(num2str(tick));
    len_pad = len_num_iters - len_tick;
    pad = '';
    for pi = 1 : len_pad
        pad = strcat('0', pad);
    end
    fig_file = strcat(path, '/', 'particle_concentration_', num2str(pt), '_', pad, num2str(tick));
    saveas(gcf, fig_file, 'pdf');
    
    % delete figure
    close(h);
    
end

end  % function save_particle_concentration

