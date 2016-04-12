%% plot certain time series and save them to disk
function [] = save_time_series(ct_range)

global num_cell_types tick population cell_type_color_map fitness_1_avg fitness_1_std fitness_n_avg fitness_n_std path;

% set up an invisible figure
h = figure('Visible', 'off');

% population by cell type (time series)
subplot(3,1,1);
hold on;
for ct = 1 : num_cell_types
    plot(1:tick, population(ct, 1:tick), 'color', cell_type_color_map(ct,:));
end
hold off;
if tick > 1
    xlim([1 tick]);
end
title('Time Evolution of Populations by Cell Type');

% local fitness by cell type (time series)
subplot(3,1,2);
plot_time_series(fitness_1_avg(ct_range,:), fitness_1_std(ct_range,:), cell_type_color_map(ct_range,:), 1, tick);
title('Time Evolution of Local Fitness by Cell Type');

% neighborhood fitness by cell type (time series)
subplot(3,1,3);
plot_time_series(fitness_n_avg(ct_range,:), fitness_n_std(ct_range,:), cell_type_color_map(ct_range,:), 1, tick);
title('Time Evolution of Neighborhood Fitness by Cell Type');

% create figure file name and save figure
fig_file = strcat(path, '/', 'time_series');
saveas(h, fig_file, 'fig');
saveas(h, fig_file, 'pdf');

% delete figure
close(h);

end  % function save_time_series

