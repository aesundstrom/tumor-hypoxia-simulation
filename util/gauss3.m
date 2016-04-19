%% create a 3D Gaussian filter using a kernel and standard deviation
%  k: kernel             ([kernel_dim_x kernel_dim_y kernel_dim_z])
%  s: standard deviation ([sigma_x sigma_y sigma_z])

function [g] = gauss3(k, s)

k = floor(k/2);
[gx, gy, gz] = ndgrid(-k(1):k(1), -k(2):k(2), -k(3):k(3));
gx = gx / s(1);
gy = gy / s(2);
gz = gz / s(3);
g = exp(-(gx .* gx + gy .* gy + gz .* gz) .* 0.5);
g = g / sum(g(:));

end  % function gauss3

