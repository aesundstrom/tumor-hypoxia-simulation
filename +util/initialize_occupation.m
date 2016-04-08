%% initialize cell types
function [] = initialize_occupation()
initialize_occupation_pattern;
occupation = permute(occupation, [2 1 3]);
end  % function initialize_occupation
