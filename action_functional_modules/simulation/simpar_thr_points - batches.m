%% Отcеивание состояний, в которые попадает система после превышения порога
% возвращаем весь вектор состояния
% Симулирование методом Монте-Карло.

%  В этой версии сохраняем данные из паралельного симулирования, а затем их
%  используем. Это значит, что с воркеров мы постоянно тягаем данные. Это
%  около 1 Гб. Т.е. лишний гигабайт туда и обратно. В следующей версии я
%  выборку сразу делаю внутри воркера, что привело к уменьшению количества
%  переылаемых данных, что в свою очередь ускорило выполнение в около 2х
%  раз.
% Вывод: старайся писать так, чтобы воркер возвращал все нужные тебе
% данные.

% Пример использования:
% -------------------------------------------------
%
% -------------------------------------------------
%
% Plot simulation results:
% -------------------------------------------------
%
%
%
% -------------------------------------------------
%
%
% Отказался от среднего между вышедшей и предыдущей точкой, потому что
% некоторые точки стали быть перед выходом. А особо разницы в (МАЛОЙ)
% выборке я не увидел
%%

function [exceed_outs, k] = simpar_thr_points(sys, time, x0, thr, ...
    exceedCount, batchSize)

exceed_outs.times = [];
exceed_outs.count = [];
exceed_outs.points = [];
exceed_outs.points_avg = [];

% double ~ 8 bytes
% 1 exceed = 1+size(sys.C,1) doubles + simcount_numbers 
% exceedCount exceeds ~ numbers of exceeds * exceedCount
% На самом деле понадобится чуть больше, потому что мы точно не знаем,
% сколько будет симуляций

% (1+size(sys.C,1))*j*8+8

% profile on

fprintf('Потребуется памяти > %d  байт\n', 8 * (2*exceedCount*size(sys.C,1)+1))

j = 0;  % фактическое количество нужных событий
k = 0;  % количество симулирвоаний

tlen = length(time);
sysBsize = size(sys.B,2);

while j < exceedCount
    batches = zeros(tlen, size(sys.C,1), batchSize);
    
%     ticBytes(gcp);
    parfor i=1:batchSize
        w = randn(tlen, sysBsize);
        out = lsim(sys, w, time, x0);
        batches(:,:,i) = out;
    end
%     tocBytes(gcp)
    
%     exceeds = getPointsAfterThr(batches, thr, time(2)-time(1));
%     exceed_outs = [exceed_outs exceeds];

%     [times, count, points, points_avg] = getPointsAfterThr(batches, thr, time(2)-time(1));
    exceeds = getPointsAfterThr(batches, thr, time(2)-time(1));
    exceed_outs.times = [exceed_outs.times exceeds.times];
    exceed_outs.count = [exceed_outs.count exceeds.count];
    exceed_outs.points = [exceed_outs.points exceeds.points];
    exceed_outs.points_avg = [exceed_outs.points_avg exceeds.points_avg];
    
    j = j + exceed_outs.count(end);
    fprintf('%d / %d \n', j, exceedCount);
    k = k + batchSize;
end

% profile viewer
end


%%
function [exceeds] = getPointsAfterThr(out, thr, dh)
% function [exceeds_times, exceeds_count, exceeds_points, exceeds_points_avg] = getPointsAfterThr(out, thr, dh)
simCount = size(out,3);

% маска выбора точек
if thr>0
    mask = squeeze( out(:,1,:) >= thr );
else
    mask = squeeze( out(:,1,:) <= thr );
end

% время первого выхода к порогу
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


