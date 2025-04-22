%% ММК по поиску траекторий
%

% clear
% clc
% close all

%%
[path2module, ~, ~] = fileparts(mfilename('fullpath'));
cd(path2module)
addpath(genpath( './../../action_functional_modules/' )) % 
clear path2module

% return

% Create a parallel pool if none exists
% if isempty(gcp)
%     parpool;
% end


%% Модель
arWindSpeed = 10;
sys = getKatWindModel(arWindSpeed, 0);

tFinal = 300;
dh = 0.01;
sim_time = 0 : dh : tFinal;

%% А-профили
thrs = linspace(2, 5, 9);


%% -----------------------------------------------------------------------
pfltime = 0:.1:15;
eps = 0.1;

pfls = getProfiles(sys, pfltime, thrs, eps);
% pfls.apathps = (pfl.apath' + ps)';  % Профиль в исходной системе координат (не в отклонениях)

fig_num = 1;
drawProfiles(pfls, sys, fig_num)

fprintf('Завершено \n')
 
if 1
    return
end

%% поиск точек пересечения порога
% моделирование Монте-Карло.

x0 = zeros(1, size(sys.A, 1));
% thr = 4.625;  % 44 сек
thr = 5;  
statevec = 1 : size(sys.A, 1);
% statevec = 1 : 6;

rng(123)

mmcCount = 1000; batchSize = 5000;
% exceedCount = 300;  batchSize = 500;  tFinal = 600;
% Всего смоделировано 3000 Elapsed time is 101.242013 seconds.

tFinal = 300;
dT = 0.01;
simTime = 0 : dT : tFinal;

% симуляция на большой длительности
tic, fprintf('Simulation... ')
% [mmc_outs, k, mmc_w] = simulate_model_parforNd_thr(sys, simTime, x0, thr, statevec, mmcCount, batchSize, @simfun);
[mmc_outs, k, mmc_w] = simulate_model_parforNd_thr(sys, simTime, x0, thr, statevec, mmcCount, batchSize);
mmcCount = size(mmc_outs,3);  % update count
toc

% траектории, которые получились
figure(5), clf, hold on, grid on, title('Траектории, которые получились в ММК (пример)')
xlabel(sys.StateName(1)), ylabel(sys.StateName(2))
plot( squeeze( mmc_outs(:,1,1:10)), squeeze( mmc_outs(:,2,1:10)) )

% %% после предыдущего симулирования у нас размерность вектора состояния 
% % не 10, а всего 6 - нет учёта внешних сил, поэтому надо пересимулировать заново
% for i=1:mmcCount
%     simOut = lsim(sys, mmc_w(:,1,i), simTime);
% end

%% выбираем участки с превышением
exdSizeTime = 15;
indent = 100;
exdTime = 0:dT:exdSizeTime-dT;

exdPathes = [];
exdPtb = [];

    [exd, timesAbove, isAbove] = findAbove_v2(size(sys.C,1)+size(sys.B,2), ...
            cat(2, mmc_outs(:,:,1), mmc_w(:,:,1)), simTime, exdSizeTime, thr, indent);
    exdPathes = cat(3, exdPathes, exd);

tic
fprintf('Выбор Участков... ')
for i=2:mmcCount
%     [exd, timesAbove, isAbove] = findAbove_v2(sys, mmc_outs(:,:,i), simTime, winSizeTime, thr, indent);
    [exd, timesAbove, isAbove] = findAbove_v2(size(sys.C,1)+size(sys.B,2), ...
            cat(2, mmc_outs(:,:,i), mmc_w(:,:,i)), simTime, exdSizeTime, thr, indent);
        
    exdPathes = cat(3, exd, exdPathes);
end

fprintf([ num2str(size(exdPathes,3)), '\n'])
toc

% clear winSizeInds 

%% График со всеми участками траекторий до выхода в КС и сравнение с А-профилем
pfl.time = 0:dT:15;
pfl = getProfile_v2(sys, pfl.time, thr);

figure(22), clf
subplot(211), hold on, grid on%, title('Инстантон в сравнении с ММК')
pl1 = line([0 exdSizeTime],[thr thr], 'Linestyle', '--', 'Color',[1 0 0]);
pl2 = plot(exdTime, exdPathes(:,1,1), 'Color',[.5 .5 .5 .1]);

for i=2:size(exdPathes,3)
    plot(exdTime, exdPathes(:,1,i), 'Color',[.5 .5 .5 .1])
end

pl3 = plot(pfl.time, pfl.apath(1,:), 'Color',[.2 .3 .8], 'Linewidth', 2);
pl4 = plot(exdTime, mean(exdPathes(:,1,:),3), 'Color',[.3 .3 .4], 'Linestyle', '-', 'Linewidth', 2);

xlabel('Время', 'FontSize', 16), ylabel(sys.StateName(1), 'FontSize', 16, 'interpreter', 'latex')
legend([pl1, pl2, pl3, pl4], 'Уровень КС', 'ММК', 'Инстантон', 'Среднее ММК', 'Location', 'SouthWest', 'FontSize', 14)
ax = gca; ax.FontSize = 14; 

subplot(212), hold on, grid on%, title('Вероятности')
% pl2 = plot(pfl.time, pfl.aaf);
pl2 = plot(pfl.time, exp(-0.001^-1*pfl.aaf), 'Linewidth', 2, 'Displayname', '$\varepsilon=0.001$');
pl2 = plot(pfl.time, exp(-0.01^-1*pfl.aaf), 'Linewidth', 2, 'Displayname', '$\varepsilon=0.01$');
pl2 = plot(pfl.time, exp(-0.1^-1*pfl.aaf), 'Linewidth', 2, 'Displayname', '$\varepsilon=0.1$');
xlabel('Время', 'FontSize', 16)
ylabel('$$P(t;\varepsilon, \varphi$$)', 'FontSize', 16, 'interpreter', 'latex')
legend('Location', 'west', 'interpreter', 'latex', 'FontSize', 14)
ax = gca; ax.FontSize = 14; 
% % График БЕЛОГО ШУМА - управления до выхода в КС и сравнение с А-профилем
% subplot(212), hold on, grid on, title('Белый шум в сравнении с ММК')
% % pl1 = line([0 winSizeTime],[thr thr], 'Linestyle', '--', 'Color',[1 0 0]);
% pl2 = plot(exdTime, exdPathes(:,end,1), 'Color',[.5 .5 .5 .1]);
% 
% for i=2:size(exdPathes,3)
%     plot(exdTime, exdPathes(:,end,i), 'Color',[.5 .5 .5 .1])
% end
% pl3 = plot(pfl.time, pfl.ptb(1,:), 'Color',[.2 .1 .7]);
% xlabel('Время', 'FontSize', 16), ylabel('\xi', 'FontSize', 16)
% legend([pl2, pl3], 'ММК', 'Инстантон', 'Location', 'South')

pathfig = ['./figs/mmc_vs_instanton_mmcCount=',num2str(mmcCount),'.fig'];
savefig(gcf, pathfig)
saveas(gcf, [pathfig,'.jpg'])


% %% График со всеми участками траекторий до выхода в КС и сравнение с А-профилем (трёхмерный)
figure(23), clf, hold on, grid on%, title('Инстантон в сравнении с ММК')
% pl1 = line3([0 winSizeTime],[thr thr], 'Linestyle', '--', 'Color',[1 0 0]);
pl2 = plot(exdPathes(:,1,1), exdPathes(:,2,1), 'Color',[.5 .5 .5 .1]);

for i=2:size(exdPathes,3)
    plot(exdPathes(:,1,i), exdPathes(:,2,i), 'Color',[.5 .5 .5 .1]);
end

plot(mean(exdPathes(:,1,:),3), mean(exdPathes(:,2,:),3), 'Color',[.3 .3 .4], 'Linestyle', '-', 'Linewidth', 2);

pl3 = plot(pfl.apath(1,:), pfl.apath(2,:), 'Color',[.2 .1 .7], 'LineWidth', 2);

xlabel(sys.StateName(1), 'FontSize', 16)
ylabel(sys.StateName(2), 'FontSize', 16)
ax = gca; ax.FontSize = 14; 
legend([pl2, pl3], 'ММК', 'Инстантон')

pathfig = ['./figs/mmc_vs_instanton_mmcCount=',num2str(mmcCount),'2d.fig'];
savefig(gcf, pathfig)
saveas(gcf, [pathfig,'.jpg'])

%%
exceeds.points = [];
exceeds.points(1,:) = squeeze(exdPathes(end,1,:));
exceeds.points(2,:) = squeeze(exdPathes(end,2,:));

figure (24), clf
hold on, grid on
xlabel(sys.StateName(1), 'interpreter', 'latex')
ylabel(sys.StateName(2), 'interpreter', 'latex')%, title('Посещаемость фазовой плоскости')
% plot(squeeze(out(:,1,:)), squeeze(out(:,2,:)), '.', 'Color', [.5 .5 .5 .2])
% plot(squeeze(out(:,1,:)), squeeze(out(:,2,:)), 'Color', [.5 .5 .5 .1])

% плотность траекторий на фазовой плоскости
x = squeeze(mmc_outs(:,1,:));
y = squeeze(mmc_outs(:,2,:));
[n,c] = hist3([x(:) y(:)], 'Nbins',[300 300]);
n = n';
cx = c{1};
cy = c{2};
% pcolor(cx,cy, n)
pcolor(cx,cy, log(n+1))
% pcolor(cx,cy, log(log(n+1)+1))
shading flat
c = colorbar;
c.Label.String = 'log(n+1)';

p2 = plot(exceeds.points(1,:),exceeds.points(2,:), '.', 'Color', [.9 .15 0], 'DisplayName', 'Точки выхода');
p1 = plot(pfls(9).apath(1,:), pfls(9).apath(2,:), 'DisplayName', 'Инстантон', 'Color', [.98 .97 .97], 'Linewidth', 2);
legend([p1, p2], 'Location', 'South', 'color', [.8 .8 .8])
ax = gca; ax.FontSize = 14; 


pathfig = ['./figs/histogram_thr=',num2str(thr),'_mmcCount=',num2str(mmcCount),'_tFinal=',num2str(tFinal),'.fig'];
savefig(gcf, pathfig)
saveas(gcf, [pathfig,'.jpg'])

%%
clc

% figHandles = findobj('Type', 'figure');
% figHandles = findall(groot, 'Type', 'figure');      % Earlier versions
figHandles = get(groot, 'Children');  % Since version R2014b

set([figHandles.CurrentAxes], 'FontSize', 14)

%% Греческие буквы через спец.символы
% pl2 = plot(pfl.time, exp(-0.001^-1*pfl.aaf), 'Displayname', [char(949) '=0.001']);
% pl2 = plot(pfl.time, exp(-0.01^-1*pfl.aaf), 'Displayname', [char(949) '=0.01']);
% pl2 = plot(pfl.time, exp(-0.1^-1*pfl.aaf), 'Displayname', [char(949) '=0.1']);
% ylabel(["$$P(t;",char(949),",$$\phi^*) \varphi$$"], 'FontSize', 16, 'interpreter', 'latex')

%% функция симулирования для ММК
function [out, w] = simfun(sys, time, x0)
w = randn(length(time), size(sys.B,2));
out = lsim(sys, w, time, x0);
end