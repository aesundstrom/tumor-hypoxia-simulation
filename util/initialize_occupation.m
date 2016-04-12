%% initialize cell types
function [] = initialize_occupation()
global x_dim y_dim z_dim mid_x mid_y mid_z mid_mid_x mid_mid_y mid_mid_z occupation;
initialize_occupation_pattern;
occupation = permute(occupation, [2 1 3]);
end  % function initialize_occupation
