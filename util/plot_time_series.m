%% plot a time series average +/- standard deviation (average curve surrounded by +/- gray patches)
%  a: average            (m time series x n ticks)
%  s: standard deviation (m time series x n ticks)
%  c: color map          (m time series x 3 {r,g,b})
%  b: begin time         (integer)
%  e: end time           (integer)

function [] = plot_time_series(a, s, c, b, e)

global tick;

domain = b : e;
gray = [0.9 0.9 0.9];
hold on;
for t = 1 : size(a,1)
    patch([domain fliplr(domain)], [a(t,domain) - s(t,domain), fliplr(a(t,domain) + s(t,domain))], gray, 'LineStyle', 'none');
    plot(domain, a(t,domain), 'color', c(t,:));
end
hold off;
freezeColors;
if tick > 1
    xlim([1 tick]);
end;

end  % function plot_time_series

