%% plot cell types and save them to disk

function [] = save_cell_types()

global occupation mid_z cell_type_color_map num_iters tick path;

% set up an invisible figure
h = figure('Visible', 'off');

% cell types (2-D slice)
imagesc(occupation(:,:,mid_z), [1 size(cell_type_color_map,1)]);
colormap(cell_type_color_map);
freezeColors;
set(gca, 'YDir', 'normal');
axis square;
title({'Cell Population' ['t=' num2str(tick)]});

% create figure file name and save figure
len_num_iters = length(num2str(num_iters));
len_tick = length(num2str(tick));
len_pad = len_num_iters - len_tick;
pad = '';
for pi = 1 : len_pad
    pad = strcat('0', pad);
end
fig_file = strcat(path, '/', 'cells_', pad, num2str(tick));
saveas(gcf, fig_file, 'pdf');

% delete figure
close(h);

end  % function save_cell_types