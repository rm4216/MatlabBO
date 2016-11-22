function [ fig ] = Plot1dim( X, fX, mean, stdDevs, xs_star, ys_star )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

fig = figure(44); clf

plot(X, mean); hold on
plot(X, mean+2*stdDevs);
hold on
plot(X, mean-2*stdDevs);
plot(X, fX, 'color', 'blue', 'linewidth', 2);
plot(xs_star, ys_star, 'o', 'color', 'red');
plot(xs_star(end, :), ys_star(end), '*', 'color', 'red');
xlim([-5, 5]);
ylim([0, 30]);
title('GP-UCB 1dim');
xlabel('x');
end

