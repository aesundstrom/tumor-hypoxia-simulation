%% create a color map using a basis and a granularity
%  c: color map basis     ([r g b])
%  n: density of gradient (integer)

function [q] = grade(c, n)

r = linspace(0,c(1),n)';
g = linspace(0,c(2),n)';
b = linspace(0,c(3),n)';
q = horzcat(r,g,b);

end  % function grade

