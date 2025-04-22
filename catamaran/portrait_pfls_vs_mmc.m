%%
% 1. Выясняем статистику движения и выхода траекторий модели из аттрактора.
% Строим гистограмму траекторий на фазовой плоскости.
% Заключение: видно, что траектории выхода "кучкуются" вдоль А-профия
% 2. Не смотрел на оценки времени -- они всёравно не похожи. Что с этим
% делать -- не ясно
% 
%%

clear
clc
% close all

%%
[path2module, ~, ~] = fileparts(mfilename('fullpath'));
cd(path2module)
addpath(genpath( './../../action_functional_modules/' )) %
clear path2module

% return

% Create a parallel pool if none exists
if isempty(gcp)
    parpool;
end


%%
arWindSpeed = 10;
sys = getKatWindModel(arWindSpeed, 0);

eps = 1;

%% Моделирование методом Монте-карло
rng(123)
simCount = 10^4;

tFinal = 300;
dh = 0.01;
sim_time = 0 : dh : tFinal;

xl0 = zeros(size(sys.A,1),1);

tic, fprintf('Simulation... ')
% [out, w] = simulate_model_parfor(sys, sim_time, xl0, simCount);
out = simulate_model_parfor2d(sys, sim_time, xl0, simCount);
% out = out + ps;  % перевод назад в координаты для нелинейной системы
toc


figure (21), clf
hold on, grid on
xlabel(sys.StateName(1)), ylabel(sys.StateName(2)), title('Посещаемость фазовой плоскости')
% plot(squeeze(out(:,1,:)), squeeze(out(:,2,:)), '.', 'Color', [.5 .5 .5 .2])
% plot(squeeze(out(:,1,:)), squeeze(out(:,2,:)), 'Color', [.5 .5 .5 .1])

% плотность траекторий на фазовой плоскости
x = squeeze(out(:,1,:));
y = squeeze(out(:,2,:));
[n,c] = hist3([x(:) y(:)], 'Nbins',[100 100]);
n = n';
cx = c{1};
cy = c{2};
% pcolor(cx,cy, n)
pcolor(cx,cy, log(n+1))
% pcolor(cx,cy, log(log(n+1)+1))
shading flat
c = colorbar;
c.Label.String = 'log(n+1)';

pathfig = ['./figs/histogram_all_eps=',num2str(eps),'_simCount=',num2str(simCount),'_tFinal=',num2str(tFinal),'.fig'];
savefig(gcf, pathfig)
saveas(gcf, [pathfig,'.jpg'])


% return


%%
pfl_time = 0 : .1 : 50;
thrs = linspace(2, 5, 9);

% thrs = [5 4 3 2 1];

thr_indx = 8;
xTf = thrs(thr_indx);

% профили на разные уровни
sys_hurw = sys;
% sys_hurw.A(end, end) = -10^-.3;
% sys_hurw.A(end, end) = -10^-5;
pfls = getProfiles(sys_hurw, pfl_time, thrs, eps);


% профили разных систем в один уровень
% sys_hurw = sys;
% p = linspace(-10^-5, -1, 10);
% sysvec = cell(1,length(p));
% for i = 1: length(p)
%     sys_hurw.A(end, end) = p(i);
%     sysvec{i} = sys_hurw;
%     norm(sysvec{i}.A,2) - norm(sys.A)
% end

% pfls = getProfiles(sysvec, pfl_time, xTf, eps);
% pfls = [];
% fprintf('Считаем А-профили... '), tic
% for i = 1:length(p)
%     pfl = getProfile_v2(sysvec{i}, sim_time, xTf);
%     pfl.pr = exp(-eps^-2 * pfl.aaf) ;
%     pfls = [pfls pfl]; 
% end
% toc




% drawPfls(sys, pfls, 10);
% drawPfls(sysvec, pfls, 10);


% % Считаем статистику
% thr_indx = 3;
% xTf = thrs(thr_indx);
exceeds = statistic(out, xTf, simCount, dh);

if exceeds.count < 1
    warning('exceeds.count < 1')
    return
end


figure(22), clf, hold on, grid on
xlabel(sys.StateName(1)), ylabel(sys.StateName(2)), title(['Плотность траекторий выхода к xTf = ', num2str(xTf)])

x = exceeds.paths.e1st;
y = exceeds.paths.e2st;
[n,c] = hist3([x(:) y(:)], 'Nbins',[150 150]);
n = n'; cx = c{1}; cy = c{2};

% pcolor(cx,cy, n)
pcolor(cx,cy, log(n+1))
% pcolor(cx,cy, log(log(n+1)))
shading flat
c = colorbar;
c.Label.String = 'log(n+1)';

p2 = plot(exceeds.points(1,:),exceeds.points(2,:), '.', 'Color', [.9 .15 0], 'DisplayName', 'Точки выхода');
% p1 = plot(pfls(thr_indx).apath(1,:), pfls(thr_indx).apath(2,:), 'DisplayName', 'А-профиль', 'Color', [.8 .1 .01], 'Linewidth', 2);

p1 = plot(pfls(8).apath(1,:), pfls(8).apath(2,:), 'DisplayName', 'Инстантон', 'Color', [.98 .97 .97], 'Linewidth', 2);
for i = [1 5 8]   % 2: length(sysvec)
    plot(pfls(i).apath(1,:), pfls(i).apath(2,:), 'DisplayName', 'Инстантон', 'Color', [.98 .97 .97], 'Linewidth', 1);
end

legend([p1, p2], 'Location', 'South', 'color', [.8 .8 .8])

pathfig = ['./figs/histogram_to_xTf=',num2str(xTf),'_simCount=',num2str(simCount),'_tFinal=',num2str(tFinal),'.fig'];
savefig(gcf, pathfig)
saveas(gcf, [pathfig,'.jpg'])


clear x y cx cy

%%
figure(24), clf
subplot(211), hold on, title(['Время выхода к \theta_{Крит}=',num2str(xTf),', N=',num2str(simCount),', T=',num2str(tFinal)])
xlabel('Время, с')
h = histogram(exceeds.times);

n = h.Values;
c = h.BinEdges;

subplot(212), hold on, grid on
% plot(pfls(thr_indx).time, pfls(thr_indx).esttime);
xlabel('Время профиля, с'), ylabel('Оценка среднего времени, с')

return


%% рисуем траектории  (требует много памяти)
figure(23), clf, hold on, grid on
xlabel(sys.StateName(1)), ylabel(sys.StateName(2)), title(['Траектории выхода к xTf = ', num2str(xTf)])
% plot(exceeds.paths.e1st, exceeds.paths.e2st, '.', 'Color', [.5 .5 .5 .2])
plot(exceeds.paths.e1st, exceeds.paths.e2st, 'Color', [.5 .5 .5 .1])
p1 = plot(pfls(thr_indx).apath(1,:), pfls(thr_indx).apath(2,:), 'DisplayName', 'А-профиль', 'Color', [.8 .1 .01], 'Linewidth', 2);
p2 = plot(exceeds.points(1,:),exceeds.points(2,:), 'o', 'Color', [.9 .5 0], 'DisplayName', 'Точки выхода');
legend([p1, p2], 'Location', 'South')

% temptime = pfls(thr_indx).time(end) - t(end);
% plot(temptime + t, exceeds.paths.e1st, 'Color', [.5 .5 .5 .1])
% p1 = plot(pfls(thr_indx).time, pfls(thr_indx).apath(1,:), 'DisplayName', 'А-профиль', 'Color', [.8 .1 .01], 'Linewidth', 2);




%% 
function [exceeds] = statistic(out, xTf, simCount, dh)
% маска выбора точек
if xTf>0
    exceeds_mask = squeeze( out(:,1,:) >= xTf );
else
    exceeds_mask = squeeze( out(:,1,:) <= xTf );
end

% время до выхода к порогу
exceeds_timeidx = zeros(1, simCount);
for i = 1:simCount
    v = find(exceeds_mask(:,i), 1);  % индексы
    if isempty(v)
        exceeds_timeidx(i) = 1;    % здесь нет ошибки, эта единица не будет учитываться
    else
        exceeds_timeidx(i) = v;
    end
end

% temp = sum(exceeds_mask);
% sum(temp>0)

exceeds_indexes = exceeds_timeidx>1;
exceeds_count = sum(exceeds_indexes);

if exceeds_count < 1
    exceeds_times = -1;
    exceeds_paths = [];
    exceeds_points = [];

else
    % финальные точки после выхода из области
    endpoint_1st = squeeze(out(:,1,exceeds_indexes));   % все индексы для выбранных  траекторий
    idx = sub2ind(size(endpoint_1st), exceeds_timeidx(exceeds_indexes), 1:size(endpoint_1st,2) );
    endpoint_2st = squeeze(out(:,2,exceeds_indexes));   % все индексы для выбранных  траекторий
    exceeds_points = [endpoint_1st(idx); endpoint_2st(idx)];   %
    
    exceeds_times = exceeds_timeidx(exceeds_indexes)*dh;
   
    exceeds_paths.e1st = endpoint_1st;
    exceeds_paths.e2st = endpoint_2st;
    
end

exceeds.times = exceeds_times ;
exceeds.count = exceeds_count ;
exceeds.points = exceeds_points;
exceeds.paths = exceeds_paths;

% Для отладки:
% figure(99);clf; plot(exceeds_endpath(1,:), exceeds_endpath(1,:), '.', 'Color', [.5 .5 .5 .2])
% figure(99);clf; hold on;plot(exceeds_endpoint(1,:), exceeds_endpoint(2,:), 'x', 'Color', [.0 .5 .05])
% figure(99);clf; hold on; grid on; plot(endpoint_1st, endpoint_2st, '.', 'Color', [.0 .5 .05])
% figure(99);clf; hold on; grid on; plot(endpoint_1st, endpoint_2st, 'Color', [.0 .5 .05 .1])
end
