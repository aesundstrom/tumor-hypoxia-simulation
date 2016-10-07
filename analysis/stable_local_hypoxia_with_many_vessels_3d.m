clear all;
clc;

global x_dim y_dim z_dim occupation seen;

% simulation constants
num_sims = 10;
x_dim    = 40;
y_dim    = 40;
z_dim    = 40;
sim_len  = 300;
sim_set  = 1 : num_sims;
edge_margin = 5;
pair_margin = 8;

% load the simulation occupation and population data

path = '~/Desktop/simulations/stable_local_hypoxia_with_many_vessels_3d';
for sim_idx = sim_set
    
    ext_path = strcat( path, '/', num2str( sim_idx-1 ) );
    
    sim_occ_file_name = strcat( ext_path, '/occupation.txt' );
    sim_occ = dlmread( sim_occ_file_name );
    occupation{sim_idx} = reshape( sim_occ, x_dim, y_dim, z_dim );
    
    sim_pop_file_name = strcat( ext_path, '/population.txt' );
    sim_pop = dlmread( sim_pop_file_name );
    population{sim_idx} = sim_pop;
    
end

% Cell types:
% 1 = vessel -- white
% 2 = empty -- blue
% 3 = alpha (viable) -- red
% 4 = beta (hypoxic) -- green
% 5 = gamma (necrotic) -- orange

%
% analyze occupation
%

for sim_idx = sim_set
    
    fprintf('Sim %d:\n', sim_idx);
    
    % find the location of the vessels
    vessels = occupation{sim_idx} == 1;
    [vessels_i, vessels_j, vessels_k] = ind2sub(size(vessels), find(vessels));

    % keep the vessels within a pair margin
    pair_keep_i{sim_idx} = [];
    pair_keep_j{sim_idx} = [];
    pair_keep_k{sim_idx} = [];
    vessel_set = 1 : numel(vessels_i);
    for ref_vessel_idx = vessel_set
        ref_i = vessels_i(ref_vessel_idx);
        ref_j = vessels_j(ref_vessel_idx);
        ref_k = vessels_k(ref_vessel_idx);
        all_pass_test = 1;
        for try_vessel_idx = vessel_set
            if try_vessel_idx == ref_vessel_idx
                continue;
            end
            try_i = vessels_i(try_vessel_idx);
            try_j = vessels_j(try_vessel_idx);
            try_k = vessels_k(try_vessel_idx);
            dist = sqrt( (ref_i - try_i)^2 + (ref_j - try_j)^2 + (ref_k - try_k)^2 );
            if dist < pair_margin
                all_pass_test = 0;
            end
        end
        if all_pass_test
            pair_keep_i{sim_idx}(numel(pair_keep_i{sim_idx})+1) = ref_i;
            pair_keep_j{sim_idx}(numel(pair_keep_j{sim_idx})+1) = ref_j;
            pair_keep_k{sim_idx}(numel(pair_keep_k{sim_idx})+1) = ref_k;
            fprintf('\t(%d,%d,%d)\n', ref_i, ref_j, ref_k);
        end
    end
    
    fprintf('\t--\n');
    
    % keep the vessels within an edge margin
    edge_keep_i{sim_idx} = [];
    edge_keep_j{sim_idx} = [];
    edge_keep_k{sim_idx} = [];
    pair_keep_set = 1 : numel(pair_keep_i{sim_idx});
    for vessel_idx = pair_keep_set
        i = pair_keep_i{sim_idx}(vessel_idx);
        j = pair_keep_j{sim_idx}(vessel_idx);
        k = pair_keep_k{sim_idx}(vessel_idx);
        if i > edge_margin && i < (y_dim - edge_margin) && j > edge_margin && j < (x_dim - edge_margin) && k > edge_margin && k < (z_dim - edge_margin)
            edge_keep_i{sim_idx}(numel(edge_keep_i{sim_idx})+1) = i;
            edge_keep_j{sim_idx}(numel(edge_keep_j{sim_idx})+1) = j;
            edge_keep_k{sim_idx}(numel(edge_keep_k{sim_idx})+1) = k;
            fprintf('\t(%d,%d,%d)\n', i, j, k);
        end
    end
    
    fprintf('\t--\n');

    % compute volume (viable and hypoxic) around each vessel
    volume_ct_34{sim_idx} = [];
    volume_ct_3{sim_idx} = [];
    edge_keep_set = 1 : numel(edge_keep_i{sim_idx});
    for vessel_idx = edge_keep_set
        i = edge_keep_i{sim_idx}(vessel_idx);
        j = edge_keep_j{sim_idx}(vessel_idx);
        k = edge_keep_k{sim_idx}(vessel_idx);
        seen{sim_idx} = zeros( x_dim, y_dim, z_dim );
        v_34 = volume_3_or_4(sim_idx, i, j, k) - 1;
        seen{sim_idx} = zeros( x_dim, y_dim, z_dim );
        v_3 = volume_3(sim_idx, i, j, k) - 1;
        volume_ct_34{sim_idx}(numel(volume_ct_34{sim_idx})+1) = v_34;
        volume_ct_3{sim_idx}(numel(volume_ct_3{sim_idx})+1) = v_3;
        fprintf('\tV_34(%d,%d,%d) = %d\n', i, j, k, v_34);
        fprintf('\tV_3(%d,%d,%d) = %d\n', i, j, k, v_3);
    end

end

% output numerical results

fprintf('\nResults:\n\n');
num_results = 0;
result_34 = [];
result_3 = [];
for sim_idx = sim_set
    edge_keep_set = 1 : numel(edge_keep_i{sim_idx});
    for vessel_idx = edge_keep_set
        i = edge_keep_i{sim_idx}(vessel_idx);
        j = edge_keep_j{sim_idx}(vessel_idx);
        k = edge_keep_k{sim_idx}(vessel_idx);
        v_34 = volume_ct_34{sim_idx}(vessel_idx);
        v_3 = volume_ct_3{sim_idx}(vessel_idx);
        fprintf('%d,%d,%d,%d,%d,%d\n', sim_idx, i, j, k, v_34, v_3);
        num_results = num_results + 1;
        result_34(num_results) = v_34;
        result_3(num_results) = v_3;
    end
end

median_v_34 = median(result_34);
median_v_3 = median(result_3);
result_v_34_mask = result_34 < 2 * median_v_34;
result_v_3_mask = result_3 < 2 * median_v_3;
mean_v_34 = mean(result_34(result_v_34_mask));
mean_v_3 = mean(result_3(result_v_3_mask));
std_v_34 = std(result_34(result_v_34_mask));
std_v_3 = std(result_3(result_v_3_mask));
cv_v_34 = std_v_34 / mean_v_34;
cv_v_3 = std_v_3 / mean_v_3;
mean_r_34 = mean((result_34(result_v_34_mask) ./ (4*pi/3)).^1/3);
mean_r_3 = mean((result_3(result_v_3_mask) ./ (4*pi/3)).^1/3);
std_r_34 = std((result_34(result_v_34_mask) ./ (4*pi/3)).^1/3);
std_r_3 = std((result_3(result_v_3_mask) ./ (4*pi/3)).^1/3);
cv_r_34 = std_r_34 / mean_r_34;
cv_r_3 = std_r_3 / mean_r_3;
fprintf('\nSummary (N=%d,%d):\n\n', num_results, sum(result_v_34_mask));
fprintf('A_34: %f, %f, %f\n', mean_v_34, std_v_34, cv_v_34);
fprintf('A_3: %f, %f, %f\n', mean_v_3, std_v_3, cv_v_3);
fprintf('r_34: %f, %f, %f\n', mean_r_34, std_r_34, cv_r_34);
fprintf('r_3: %f, %f, %f\n', mean_r_3, std_r_3, cv_r_3);

%
% analyze population
%

% plot the indivual simulation time series for cell type 3 and 4 populations
figure;
hold on;
for sim_idx = sim_set
    for pop_idx = [3 4 5]
        plot( 1:sim_len, population{sim_idx}(pop_idx,:)', 'Color', [0.8 0.8 0.8] );
    end
end

% define orange color
orange = [1 0.5 0];

% compute and plot the mean time series for each cell type
sum_pop_3 = zeros(1, sim_len);
sum_pop_4 = zeros(1, sim_len);
sum_pop_5 = zeros(1, sim_len);
for sim_idx = sim_set
    sum_pop_3 = sum_pop_3 + population{sim_idx}(3,:);
    sum_pop_4 = sum_pop_4 + population{sim_idx}(4,:);
    sum_pop_5 = sum_pop_5 + population{sim_idx}(5,:);

end
avg_pop_3 = sum_pop_3 / num_sims;
avg_pop_4 = sum_pop_4 / num_sims;
avg_pop_5 = sum_pop_5 / num_sims;
plot( avg_pop_3', 'r', 'LineWidth', 4);
plot( avg_pop_4', 'g', 'LineWidth', 4);
plot( avg_pop_5', 'Color', orange, 'LineWidth', 4);
hold off;
title( 'Population Size Versus Time over 10 Simulations' );
xlabel( 'time' );
ylabel( 'number of cells' );

% compute the standard deviation in the time series
ssd_pop_3 = zeros(1, sim_len);
ssd_pop_4 = zeros(1, sim_len);
ssd_pop_5 = zeros(1, sim_len);
for sim_idx = sim_set
    ssd_pop_3 = ssd_pop_3 + (population{sim_idx}(3,:) - avg_pop_3).^2;
    ssd_pop_4 = ssd_pop_4 + (population{sim_idx}(4,:) - avg_pop_4).^2;
    ssd_pop_5 = ssd_pop_5 + (population{sim_idx}(5,:) - avg_pop_5).^2;
end
std_pop_3 = sqrt( ssd_pop_3 / num_sims );
std_pop_4 = sqrt( ssd_pop_4 / num_sims );
std_pop_5 = sqrt( ssd_pop_5 / num_sims );

% plot the standard deviation in the time series
figure;
hold on;
plot( std_pop_3', 'r' );
plot( std_pop_4', 'g' );
plot( std_pop_5', 'Color', orange );
hold off;
title( 'Standard Deviation in Population Size Versus Time over 10 Simulations' );
xlabel( 'time' );
ylabel( 'number of cells' );

% compute the coefficient of variation in the time series
var_pop_3 = std_pop_3 ./ avg_pop_3;
var_pop_4 = std_pop_4 ./ avg_pop_4;
var_pop_5 = std_pop_5 ./ avg_pop_5;

% plot the coefficient of variation in the time series
figure;
hold on;
plot( var_pop_3', 'r' );
plot( var_pop_4', 'g' );
plot( var_pop_5', 'Color', orange );
hold off;
title( 'Coefficient of Variation in Population Size Versus Time over 10 Simulations' );
xlabel( 'time' );
ylabel( 'CV' );