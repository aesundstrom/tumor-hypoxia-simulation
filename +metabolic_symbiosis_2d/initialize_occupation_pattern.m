%% initialize occupation pattern
function [] = initialize_occupation_pattern()
occupation(:, 1,             mid_z) = 1;  % place a layer of vessel at the bottom
occupation(:, 2:mid_y,       mid_z) = 4;  % place a stack of layers of betas just above the vessels
occupation(:, mid_y+1:y_dim, mid_z) = 3;  % place a stack of layers of alphas just above the alphas
end  % function initialize_occupation_pattern
