%% Отcеивание состояний, в которые попадает система после превышения порога
% возвращаем весь вектор состояния
% Симулирование методом Монте-Карло.
% Этот метод учитывает только первые моменты выхода.

% Вывод: старайся писать так, чтобы воркер возвращал все нужные тебе
% данные, вместо того, чтобы копить в нём.

% Пример использования:
% -------------------------------------------------
%
% -------------------------------------------------
%
% Plot simulation results:
% -------------------------------------------------
%
% -------------------------------------------------
%
%
% Почти Отказался от среднего между вышедшей и предыдущей точкой, потому что
% некоторые точки стали быть перед выходом. А особо разницы в (МАЛОЙ)
% выборке я не увидел. В большой выборке разница есть.

% Катамаран, 300шт на thr=5, Tf=600
% Elapsed time is 1117.625550 seconds. ~ 18.6 min
% 
% Катамаран, 1000шт на thr=4.5, Tf=600/ (Всего смоделировано траекторий 67000 )
% Elapsed time is 467.670598 seconds.
% 
% Хищник-жертва, 1000, thr=0.6, Tf=600 (total 33000)
% Elapsed time is 262.183232 seconds.
% 
% Хищник-жертва, 1000, thr=0.7, Tf=600 (total 190000)
% Elapsed time is 26 minutes.
% 
%%

function [exceed_outs, simCount] = simpar_thr_points(sys, time, x0, thr, ...
    exceedCount, batchSize)
fprintf('Потребуется памяти: \n')
fprintf('\tдля результата > %d  байт\n', 8 * (2*exceedCount*size(sys.C,1)+exceedCount))
% fprintf('\tдля моделирования > %d  байт\n', 8 * (2*exceedCount*size(sys.C,1)+exceedCount))

times = [];
points = [];
points_avg = [];

j = 0;  % фактическое количество нужных событий
simCount = 0;  % количество симулирвоаний

tlen = length(time);
sysBsize = size(sys.B,2);

% profile on

while j < exceedCount
%     ticBytes(gcp);

	tic
    parfor i=1:batchSize
        w = randn(tlen, sysBsize);
        out = lsim(sys, w, time, x0);
    
        exceeds = getPointsAfterThr(out, thr, time(2)-time(1));
        if exceeds.count > 0
            times = [times exceeds.times];
            points = [points exceeds.points];   % append column
            points_avg = [points_avg exceeds.points_avg];   % append column
            j = j + 1;
        end
    end
    bTime = toc;
%     bTime = (bTime+toc)/j;
%     tocBytes(gcp)
    
    fprintf('%d/%d; \t batchTime = %.2f seconds \n', j, exceedCount, bTime);
    simCount = simCount + batchSize;
end

% profile viewer

exceed_outs.times = times;
exceed_outs.points = points;
exceed_outs.points_avg = points_avg;
end


%%
function [exceeds] = getPointsAfterThr(out, thr, dh)
simCount = size(out,3);

% маска выбора точек
if thr>0
    mask = squeeze( out(:,1,:) >= thr );
else
    mask = squeeze( out(:,1,:) <= thr );
end

% время первого выхода к порогу (!!!!!!!)
timeidx = zeros(1, simCount);
for i = 1:simCount
    v = find(mask(:,i), 1);  % индексы
    if isempty(v)
        timeidx(i) = 1;    % здесь нет ошибки, эта единица не будет учитываться
    else
        timeidx(i) = v;
    end
end

indexes = timeidx>1;
count = sum(indexes);

if count < 1
    times = -1;
    points = [];
    points_avg = [];
else
    % все индексы для выбранных  траекторий
    endpoint_1st = squeeze(out(:,1,indexes));
    idx = sub2ind(size(endpoint_1st), timeidx(indexes), 1:size(endpoint_1st,2) );
    
    points = [];
    points_avg = [];
    % остальные размерности вектора состяония
    for i=1: size(out,2)
        endpoint = squeeze(out(:,i,indexes));

        points = [points; endpoint(idx)];
        points_avg = [points_avg; (endpoint(idx) + endpoint(idx-1))/2];
    end
    
    times = timeidx(indexes)*dh;
end

exceeds.times = times ;
exceeds.count = count ;
exceeds.points = points;
exceeds.points_avg = points_avg;

% Для отладки:
% figure(99),clf, hold on, grid on
% plot(endpoint_1st(idx), endpoint_2st(idx), '.', 'Displayname', 'after thr')
% plot(endpoint_1st(idx-1), endpoint_2st(idx-1), 'o', 'Displayname', 'before thr')
% plot(exceeds_points(1,:), exceeds_points(2,:), 'x', 'Displayname', 'average')
% line([thr thr],[min( exceeds_points(2,:)) max( exceeds_points(2,:))], 'Displayname', 'thr')
% legend
end


