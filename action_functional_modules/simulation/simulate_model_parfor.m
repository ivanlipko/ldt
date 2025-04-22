%% Симулирование методом Монте-Карло
% 
% Пример использования: 
% -------------------------------------------------
% x0 = [0 0]';   % начальное состояние
% time = 0 : 0.01 : 3;
% simCount = 10000;
% 
% [path2module, ~, ~] = fileparts(mfilename('fullpath'));
% cd(path2module)
% addpath(genpath( '../action_functional_modules' ))
% sys = getOscModel();
% 
% tic, fprintf('Simulation... ')
% [outs, w] = simulate_model(sys, x0, simCount, time, 3)
% toc
% -------------------------------------------------
% 
% Plot simulation results:
% -------------------------------------------------
% thr = .5*10^-5;  % уровень превышения.
% strLabel.Time = 'Время, с';
% tic, fprintf('Plot...')
% figure(1), clf, hold on, grid on
% plot(time, squeeze(outs(:,1,:)), 'Color', [.5 .5 .5 .25])
% title('Все траектории')
% line([T0 Tf],[thr thr], 'Color', [.8 .05 .1])
% xlabel(strLabel.Time)
% toc
% -------------------------------------------------
%
% save 'simulation.mat' outs w time
%
%
% Симулирование с помощью белого шума
% 
% https://www.gaussianwaves.com/2013/11/simulation-and-analysis-of-white-noise-in-matlab/
%
%%

function [outs, w] = simulate_model_parfor(sys, time, x0, simCount)
w = randn(length(time), size(sys.B,2), simCount);
% outs = zeros(length(time), size(sys.C,2), simCount);
outs = zeros(length(time), size(sys.C,1), simCount);

parfor i=1:simCount
    outs(:,:,i) = lsim(sys, w(:,:,i), time, x0);
end

end
