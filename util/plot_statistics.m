%% plot the statistics of interest

function [] = plot_statistics()

global num_cell_types num_particle_types occupation cell_type_color_map x_dim y_dim z_dim mid_z diameter_x diameter_y diameter_z tick wallclock plot_3d fitness_1 fitness_1_avg fitness_1_std fitness_n fitness_n_avg fitness_n_std population epc concentration output_results output_cell_types output_time_series output_local_fitness_3d output_local_fitness output_neighborhood_fitness output_particle_concentration;

% define some useful constants
gray_256    = grade([1 1 1], 256);
num_rows    = 6;
num_ct_cols = (num_cell_types - 2) + 1;  % no columns for vessel and empty cell types + 1 column for corresponding time series
num_cols    = max(num_ct_cols, num_particle_types);
del_cols    = num_cols - num_ct_cols;
ct_range    = 3 : num_cell_types;
pt_range    = 1 : num_particle_types;
cursor      = 1;

% set up the figure
figure(1);
clf;

% cell types (2-D slice)
subplot(num_rows, num_cols, cursor);
imagesc(occupation(:,:,mid_z), [1 size(cell_type_color_map,1)]);
colormap(cell_type_color_map);
freezeColors;
set(gca, 'YDir', 'normal');
axis square;
cursor = cursor + 1;

% performance (time series)
subplot(num_rows, num_cols, cursor);
domain = 1:tick-1;
wc     = wallclock(domain);
wc_avg = mean(wc);
wc_std = std(wc);
plot(domain, wc);
freezeColors;
title(['avg=' sprintf('%3.3f', wc_avg) ', std=' sprintf('%3.3f', wc_std)]);
axis square;
cursor = cursor + 1;

% pad, if necessary
pad_cursor();

% cell type 3D plot
for ct = ct_range
    if plot_3d
        subplot(num_rows, num_cols, cursor);
        occ_ct = occupation == ct;
        fit_ct = fitness_1(occ_ct);
        min_fit_ct = min(fit_ct);
        max_fit_ct = max(fit_ct);
        if min_fit_ct < max_fit_ct
            norm_fit_ct = (fit_ct - min_fit_ct) ./ (max_fit_ct - min_fit_ct);
        else
            norm_fit_ct = 0;
        end
        color_density = 256;
        color_grade = grade(cell_type_color_map(ct,:), color_density);
        color_idx = floor((color_density - 1) * norm_fit_ct) + 1;
        scatter_colors = color_grade(color_idx,:);
        [x,y,z] = ind2sub(size(occ_ct), find(occ_ct));
        % note the axis order for plotting: <y,x,z>
        scatter3(y, x, z, 20, scatter_colors, 'filled', 's');
        colormap(grade(cell_type_color_map(ct,:), color_density));
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
    end
    cursor = cursor + 1;
end

% population by cell type (time series)
subplot(num_rows, num_cols, cursor);
hold on;
for ct = 1 : num_cell_types
    plot(1:tick, population(ct, 1:tick), 'color', cell_type_color_map(ct,:));
end
hold off;
freezeColors;
if tick > 1
    xlim([1 tick]);
end
axis square;
cursor = cursor + 1;

% pad, if necessary
pad_cursor();

% cell type local fitness (2-D slice)
for ct = ct_range
    subplot(num_rows, num_cols, cursor);
    imagesc((occupation(:,:,mid_z) == ct) .* fitness_1(:,:,mid_z));
    colormap(grade(cell_type_color_map(ct,:), 256));
    cbfreeze(colorbar);
    freezeColors;
    set(gca, 'YDir', 'normal');
    axis square;
    cursor = cursor + 1;
end

% local fitness by cell type (time series)
subplot(num_rows, num_cols, cursor);
plot_time_series(fitness_1_avg(ct_range,:), fitness_1_std(ct_range,:), cell_type_color_map(ct_range,:), 1, tick);
axis square;
cursor = cursor + 1;

% pad, if necessary
pad_cursor();

% cell type neighborhood fitness (2-D slice)
for ct = ct_range
    subplot(num_rows, num_cols, cursor);
    imagesc((occupation(:,:,mid_z) > 1) .* fitness_n(:,:,mid_z,ct));
    colormap(grade(cell_type_color_map(ct,:), 256));
    cbfreeze(colorbar);
    freezeColors;
    set(gca, 'YDir', 'normal');
    axis square;
    cursor = cursor + 1;
end

% neighborhood fitness by cell type (time series)
subplot(num_rows, num_cols, cursor);
plot_time_series(fitness_n_avg(ct_range,:), fitness_n_std(ct_range,:), cell_type_color_map(ct_range,:), 1, tick);
axis square;
cursor = cursor + 1;

% pad, if necessary
pad_cursor();

% cell type x,y,z-extent
for ct = ct_range
    subplot(num_rows, num_cols, cursor);
    hold on;
    % note the axis order for plotting: <y,x,z>
    plot(1:tick, diameter_y(ct, 1:tick), 'r');
    plot(1:tick, diameter_x(ct, 1:tick), 'g');
    plot(1:tick, diameter_z(ct, 1:tick), 'b');
    hold off;
    if tick > 1
        xlim([1 tick]);
    end
    axis square;
    cursor = cursor + 1;
end

% Euler-Poincare characteristic by cell type (time series)
subplot(num_rows, num_cols, cursor);
hold on;
for ct = 1 : num_cell_types
    plot(1:tick, epc(ct, 1:tick), 'color', cell_type_color_map(ct,:));
end
hold off;
freezeColors;
axis square;
cursor = cursor + 1;

% pad, if necessary
pad_cursor();

% particle type concentrations (2-D slice)
for pt = pt_range
    subplot(num_rows, num_cols, cursor);
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
    title(['min' min_rel sprintf('%3.3f', min_con_pt) ', max' max_rel sprintf('%3.3f', max_con_pt)]);
    axis square;
    xlim([1 x_dim]);
    ylim([1 y_dim]);
    cursor = cursor + 1;
end

% if outputting results, then create and save a stand-alone figure
if output_results
    if output_cell_types
        save_cell_types();
    end
    if output_time_series
        save_time_series(ct_range);
    end
    if output_local_fitness_3d
        save_local_fitness_3d(ct_range);
    end
    if output_local_fitness
        save_local_fitness(ct_range);
    end
    if output_neighborhood_fitness
        save_neighborhood_fitness(ct_range);
    end
    if output_particle_concentration
        save_particle_concentration(gray_256, pt_range)
    end
end

    function [rc] = row_cursor(cursor, num_cols)
        rc = mod(cursor, num_cols);
        if rc == 0
            rc = num_cols;
        end
    end  % function row_cursor

    function [] = pad_cursor()
        for p = 1 : (num_cols - row_cursor(cursor - 1, num_cols))
            cursor = cursor + 1;
        end
    end  % function pad_cursor

end  % function plot_statistics
