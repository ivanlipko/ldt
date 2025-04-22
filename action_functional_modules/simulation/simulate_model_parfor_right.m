%% Симулирование методом Монте-Карло
%
% -------------------------------------------------
% x0 = [10 0]';   % начальное состояние
% time = 0 : 0.01 : 3;
% simCount = 10;
% rng(124)
% 
% [path2module, ~, ~] = fileparts(mfilename('fullpath'));
% cd(path2module)
% addpath(genpath( '../action_functional_modules' ))
% sys = getOscModel();
% 
% tic, fprintf('Simulation... ')
% ptb = 10*randn(length(time), size(sys.B,2), simCount);
% outs = simulate_model_parfor_right             (sys, sys.B, x0, time, ptb, @rightdpt, simCount);
% toc
% 
% % Plot simulation results:
% % -------------------------------------------------
% tic, fprintf('Plot...')
% figure(1), clf, hold on, grid on
% plot(time, squeeze(outs(:,1,:)), 'Color', [.5 .5 .5 .25])
% title('Все траектории')
% xlabel('Время, с')
% toc
% 
% % -------------------------------------------------
% function dx = rightdpt(t, x, A,ptb)
% % dx = A * x + B*cntrl + G*ptb;
% 
% dx = A * x + ptb;
% end
% 
% -------------------------------------------------
%
% save 'simulation.mat' outs w time
%
%
% Симулирование с помощью белого шума
%
% https://www.gaussianwaves.com/2013/11/simulation-and-analysis-of-white-noise-in-matlab/
%
%
%%
% inSignal = Time X State X SimCount

function [outs] = simulate_model_parfor_right(nsys, xstate0, time, inSignal, right, simCount)
if nargin == 7
    exit('simulate_model_parfor_right изменён!  Используй simulate_model_parfor_rightsys')
end

dt = time(2) - time(1);
outs = zeros(length(time), nsys, simCount);

parfor c = 1:simCount
    state = xstate0;
    i = 1;    
    out = zeros(length(time), nsys);  %  инициализация матрицы
    for t = time
        out(i,:) = state;    % присвоение н.у.
        state = rkStep(t, dt, state, right, inSignal(i,:,c));
        i = i + 1;
    end
    outs(:,:,c) = out;
end

end



%% Один шаг метода Рунге-Кутта 4-го порядка
function f = rkStep(t, h, X, right, inSignal)

h2 = 0.5*h;
h6 = 0.166666666*h;
Fs = right(t,X, inSignal);
t = t + h2;

Xr = X + h2 * Fs;

F = right(t,Xr, inSignal);

s=F;
Fs=Fs+s+s;
Xr=X+h2*s;

F = right(t,Xr, inSignal);
t=t+h2;

s=F;
Fs=Fs+s+s;
Xr=X+h*s;

F = right(t,Xr, inSignal);

X=X+h6*(Fs+F);

f=X;
end
