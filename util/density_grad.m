%% create a random 3D binary matrix using a gradient function for density along any combination of axis directions
%  xd   : x dimension
%  yd   : y dimension
%  zd   : z dimension
%  f    : gradient function (f(upper)...f(lower) guaranteed to be in the range [0,1])
%  d_beg: domain lower bound
%  d_end: domain upper bound
%  dir  : direction of gradient (3-element binary vector)
function [dm] = density_grad(xd, yd, zd, f, d_beg, d_end, dir)
dm = [];

domain = linspace(d_beg,d_end,xd);
target = f(domain);
dx = [];
for z_idx = 1 : zd
    zm = [];
    for x_idx = 1 : xd
        zm = horzcat(zm, rand(yd,1) < target(x_idx));
    end
    dx(:,:,z_idx) = zm;
end

domain = linspace(d_beg,d_end,yd);
target = f(domain);
dy = [];
for z_idx = 1 : zd
    zm = [];
    for y_idx = 1 : yd
        zm = vertcat(zm, rand(1,xd) < target(y_idx));
    end
    dy(:,:,z_idx) = zm;
end

domain = linspace(d_beg,d_end,zd);
target = f(domain);
dz = [];
for z_idx = 1 : zd
    zm = rand(xd,yd) < target(z_idx);
    dz(:,:,z_idx) = zm;
end

if     isequal(dir, [0 0 0])
    dm = rand(xd,yd,zd);
elseif isequal(dir, [0 0 1])
    dm = dz;
elseif isequal(dir, [0 1 0])
    dm = dy;
elseif isequal(dir, [0 1 1])
    dm = dy & dz;
elseif isequal(dir, [1 0 0])
    dm = dx;
elseif isequal(dir, [1 0 1])
    dm = dx & dz;
elseif isequal(dir, [1 1 0])
    dm = dx & dy;
elseif isequal(dir, [1 1 1])
    dm = dx & dy & dz;
end
end  % function density_grad

