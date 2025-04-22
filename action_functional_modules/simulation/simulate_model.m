% Симулирование с помощью белого шума
% 
% https://www.gaussianwaves.com/2013/11/simulation-and-analysis-of-white-noise-in-matlab/
%


function [sim, w] = simulate_model(sys, time, x0)
% w = wgn(size(time,2), size(sys.B,2), 0);
w = randn(size(time,2), size(sys.B,2));
sim = lsim(sys, w, time, x0);
