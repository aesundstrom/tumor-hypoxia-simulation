%% initialize occupation pattern
occupation(:,     :,     :    ) = 2;  % initialize with empty cells
occupation(mid_x, mid_y, mid_z) = 3;  % place an alpha smack in the middle
for i = 1 : 100  % randomly place 100 vessels
    rx = mid_x;
    while rx == mid_x
        rx = floor(rand * x_dim) + 1;
    end
    ry = mid_y;
    while ry == mid_y
        ry = floor(rand * y_dim) + 1;
    end
    rz = mid_z;
    while rz == mid_z
        rz = floor(rand * z_dim) + 1;
    end
    occupation(rx, ry, rz) = 1;
end
