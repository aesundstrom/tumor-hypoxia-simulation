clear all;
clc;

% simulation constants
num_sims = 10;
x_dim    = 40;
y_dim    = 40;
sim_len  = 100;
sim_set  = 1 : num_sims;

% load the simulation occupation and population data

path = '../simulations/metabolic_symbiosis_2d';
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
% 1 = vessel -- while
% 2 = empty -- blue
% 3 = alpha (hypoxia) -- red
% 4 = beta (aerobic) -- green

%
% analyze occupation
%

% compute spatial frequency of each occupation map

sum_mag_occ = zeros(y_dim-1, x_dim);
for sim_idx = sim_set
    
    % strip off the top row of vessels
    rel_occ = occupation{sim_idx}(2:y_dim, :);
    
    % create a binary occupation map
    bin_occ = rel_occ == 3;
    
    % compute the fft of the binary occupation map
    fft_occ = fft2( bin_occ );
    
    % compute the magnitude of the fft
    mag_occ{sim_idx} = abs( fftshift( fft_occ ) );

    % add the magnitude of the fft with previous sims, for averaging later
    sum_mag_occ = sum_mag_occ + mag_occ{sim_idx};
    
end

% compute the mean magnitude of the fft
avg_mag_occ = sum_mag_occ / num_sims;

% show it, scaling up to the second largest numerical value (the first is always the center position in the frequency portrait, which throws off the scaling)
srt_mag_occ = sort( reshape( avg_mag_occ, 1, x_dim * (y_dim-1) ), 'descend' );
figure; imshow(avg_mag_occ, [0 srt_mag_occ(2)]); colormap gray; colorbar;
title( strcat( 'Mean FFT2 Magnitude over 10 Simulations' ) );

% sum the squared differences from the mean magnitude of the fft, for computing the stardard devation later
ssd_mag_occ = zeros(y_dim-1, x_dim);
for sim_idx = sim_set
    ssd_mag_occ = ssd_mag_occ + (mag_occ{sim_idx} - avg_mag_occ).^2;
end

% compute and show the standard deviation in the magnitude of the fft
std_mag_occ = sqrt( ssd_mag_occ / num_sims );
figure; imshow(std_mag_occ, []); colormap gray; colorbar;
title( strcat( 'Standard Deviation in FFT2 Magnitude over 10 Simulations' ) );

% compute and show the coefficient of variation in the magnitude of the fft
var_mag_occ = std_mag_occ ./ avg_mag_occ;
figure; imshow(var_mag_occ, []); colormap gray; colorbar;
title( strcat( 'Coefficient of Variation in FFT2 Magnitude over 10 Simulations' ) );

%
% analyze population
%

% plot the indivual simulation time series for cell type 3 and 4 populations
figure;
hold on;
for sim_idx = sim_set
    for pop_idx = [3 4]
        plot( 1:sim_len, population{sim_idx}(pop_idx,:)', 'Color', [0.8 0.8 0.8] );
    end
end

% compute and plot the mean time series for each cell type
sum_pop_3 = zeros(1, sim_len);
sum_pop_4 = zeros(1, sim_len);
for sim_idx = sim_set
    sum_pop_3 = sum_pop_3 + population{sim_idx}(3,:);
    sum_pop_4 = sum_pop_4 + population{sim_idx}(4,:);
end
avg_pop_3 = sum_pop_3 / num_sims;
avg_pop_4 = sum_pop_4 / num_sims;
plot( avg_pop_3', 'r', 'LineWidth', 4);
plot( avg_pop_4', 'g', 'LineWidth', 4);
hold off;
title( 'Population Size Versus Time over 10 Simulations' );
xlabel( 'time' );
ylabel( 'number of cells' );

% compute the standard deviation in the time series
ssd_pop_3 = zeros(1, sim_len);
ssd_pop_4 = zeros(1, sim_len);
for sim_idx = sim_set
    ssd_pop_3 = ssd_pop_3 + (population{sim_idx}(3,:) - avg_pop_3).^2;
    ssd_pop_4 = ssd_pop_4 + (population{sim_idx}(4,:) - avg_pop_4).^2;
end
std_pop_3 = sqrt( ssd_pop_3 / num_sims );
std_pop_4 = sqrt( ssd_pop_4 / num_sims );

% plot the standard deviation in the time series
figure;
hold on;
plot( std_pop_3', 'r' );
plot( std_pop_4', 'g' );
hold off;
title( 'Standard Deviation in Population Size Versus Time over 10 Simulations' );
xlabel( 'time' );
ylabel( 'number of cells' );

% compute the coefficient of variation in the time series
var_pop_3 = std_pop_3 ./ avg_pop_3;
var_pop_4 = std_pop_4 ./ avg_pop_4;

% plot the coefficient of variation in the time series
figure;
hold on;
plot( var_pop_3', 'r' );
plot( var_pop_4', 'g' );
hold off;
title( 'Coefficient of Variation in Population Size Versus Time over 10 Simulations' );
xlabel( 'time' );
ylabel( 'CV' );
