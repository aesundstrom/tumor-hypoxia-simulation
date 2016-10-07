function [r] = area_3_or_4( s, i, j )

global x_dim y_dim occupation seen;

if i < 1 || i > y_dim || j < 1 || j > x_dim
    r = 0;
else
    ct = occupation{s}(i,j);
    if seen{s}(i,j)
        r = 0;
    elseif ct == 5
        r = 0;
    elseif ct == 1 || ct == 3 || ct == 4
        seen{s}(i,j) = 1;
        r = 1 + area_3_or_4(s, i-1, j-1) ...
              + area_3_or_4(s, i-1, j  ) ...
              + area_3_or_4(s, i-1, j+1) ...
              + area_3_or_4(s, i  , j-1) ...
              + area_3_or_4(s, i  , j+1) ...
              + area_3_or_4(s, i+1, j-1) ...
              + area_3_or_4(s, i+1, j  ) ...
              + area_3_or_4(s, i+1, j+1) ;
    else
        r = 0;
    end
end

return;

end