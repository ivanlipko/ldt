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
% outs = simulate_model_parfor_right(sys, sys.B, x0, time, ptb, @rightdpt, simCount);
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

function [outs] = simulate_model_parfor_right(sys, sysG, x0, time, ptb, right, count)
dt = time(2) - time(1);

outs = zeros(length(time), size(sys.C,2), count);

sysA = sys.A;
% sysB = sys.B;

parfor c = 1:count
    state = x0;
    i = 1;
    
    out = zeros(length(time), size(sys.C,2));
    for t = time
        out(i,:) = state;
        state = rkStep(t, dt, state, right, sysA, sysG * ptb(i,:,c));
%         state = rkStep(t, dt, state, right, sysA, sysG * ptb(i,:,c))';

        i = i + 1;
    end
    outs(:,:,c) = out;
end

end



%% Один шаг метода Рунге-Кутта 4-го порядка
function f = rkStep(t, h, X, right, A, ptb)

h2 = 0.5*h;
h6 = 0.166666666*h;
Fs = right(t,X, A,ptb);
t = t + h2;

Xr = X + h2 * Fs;

F = right(t,Xr, A,ptb);

s=F;
Fs=Fs+s+s;
Xr=X+h2*s;

F = right(t,Xr, A,ptb);
t=t+h2;

s=F;
Fs=Fs+s+s;
Xr=X+h*s;

F = right(t,Xr, A,ptb);

X=X+h6*(Fs+F);

f=X;
end
