clear all;
clc;

% simulation constants
num_sims = 10;
x_dim    = 40;
y_dim    = 40;
sim_len  = 300;
sim_set  = 1 : num_sims;

% load the simulation occupation and population data

path = '~/Desktop/simulations/tumor_stroma_signaling_2d';
for sim_idx = sim_set
    
    ext_path = strcat( path, '/', num2str( sim_idx-1 ) );
    
    sim_occ_file_name = strcat( ext_path, '/occupation.txt' );
    sim_occ = dlmread( sim_occ_file_name );
    occupation{sim_idx} = reshape( sim_occ, x_dim, y_dim );
    
    sim_pop_file_name = strcat( ext_path, '/population.txt' );
    sim_pop = dlmread( sim_pop_file_name );
    population{sim_idx} = sim_pop;
    
end

% Cell types:
% 1 = vessel -- white (N/A)
% 2 = empty -- blue
% 3 = alpha (epithelial) -- red
% 4 = beta (fibroblast) -- green
% 5 = gamma (tumor) -- orange
% 6 = delta (inert) -- purple

%
% analyze population
%

% plot the indivual simulation time series for cell type 5 population
figure;
hold on;
for sim_idx = sim_set
    for pop_idx = [5]
        plot( 1:sim_len, population{sim_idx}(pop_idx,:)', 'Color', [0.8 0.8 0.8] );
    end
end

% define orange color
orange = [1 0.5 0];

% compute and plot the mean time series for each cell type
sum_pop_5 = zeros(1, sim_len);
for sim_idx = sim_set
    sum_pop_5 = sum_pop_5 + population{sim_idx}(5,:);
end
avg_pop_5 = sum_pop_5 / num_sims;
plot( avg_pop_5', 'Color', orange, 'LineWidth', 4);
hold off;
title( 'Population Size Versus Time over 10 Simulations' );
xlabel( 'time' );
ylabel( 'number of cells' );

% compute the standard deviation in the time series
ssd_pop_5 = zeros(1, sim_len);
for sim_idx = sim_set
    ssd_pop_5 = ssd_pop_5 + (population{sim_idx}(5,:) - avg_pop_5).^2;
end
std_pop_5 = sqrt( ssd_pop_5 / num_sims );

% plot the standard deviation in the time series
figure;
hold on;
plot( std_pop_5', 'Color', orange );
hold off;
title( 'Standard Deviation in Population Size Versus Time over 10 Simulations' );
xlabel( 'time' );
ylabel( 'number of cells' );

% compute the coefficient of variation in the time series
var_pop_5 = std_pop_5 ./ avg_pop_5;

% plot the coefficient of variation in the time series
figure;
hold on;
plot( var_pop_5', 'Color', orange );
hold off;
title( 'Coefficient of Variation in Population Size Versus Time over 10 Simulations' );
xlabel( 'time' );
ylabel( 'CV' );
