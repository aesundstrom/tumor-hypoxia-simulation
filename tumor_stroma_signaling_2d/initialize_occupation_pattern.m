%% initialize occupation pattern
rg = density_grad(x_dim, mid_mid_y-1, 1, @(x) x.^7, 0.0, 0.8, [1 0 0]);
rg(rg == 1) = 4;
rg(rg == 0) = 6;
occupation(:,        :,                      :    ) = 2;    % initialize with empty cells
occupation(1:x_dim,  mid_mid_y-1:mid_mid_y,  mid_z) = 3;    % lay down a layer of alpha cells
occupation(mid_x,    mid_mid_y,              mid_z) = 5;    % place a gamma cell in the middle
occupation(1:x_dim,  1:mid_mid_y-1,          mid_z) = rg';  % lay down a deeper sub-epithelial region of progressively dense beta cells mixed with delta cells
