% Симулирование модели методом Монте-Карло
%
% оценка вероятности времени первого выхода из области
% P (t < T), t = {t: x > XThr}
% 
%  (!!!) надо быть аккуратным, когда меняешь границы времени. Оно хорошо
%  работает, когда такт дискретности, длительность совпадает в Методе
%  Монте-Карло и в методе стрельбы

[path2module, ~, ~] = fileparts(mfilename('fullpath'));
cd(path2module)
addpath(genpath( '../../action_functional_modules' ))

%% Симулирование. Расчёт
time = 0 : 0.01 : 1;
simCount = 5000;

x0 = [0];
sys = get1Dsys();

% [outs, w] = simulate_model(sys, x0, simCount, time, 3);
[outs, w] = simulate_model_parfor(sys, time, x0, simCount);

% outDet = lsim(sys, zeros(size(time,2),1), time, x0);   % невозмущённое движение

outs = squeeze(outs);
w = squeeze(w);


%% ФД S = int(dot_x - b(x)) и вероятности. Расчёт
tic, fprintf('ФД... ')

afsMC =  zeros(simCount,1);
parfor i = 1:simCount
    afsMC(i) = getAF1D(sys, outs(:,i));
end
toc

% probs = exp(- afsMC / sys.B^2);


%% ТРаектории. Рисунки
tic, fprintf('Plot...')

figure(1), clf, hold on, grid on
plot(time, outs, 'Color', [.5 .5 .5 .1])
title('Траектории'), xlabel('Время, с')

toc

%% Зависимость положения концов траекторий и ФД методом Монте-Карло
figure(4), clf, hold on, grid on
plot(outs(end,:), afsMC, '.', 'Color', [.5 .5 .5 .5])
% plot(outDet(end), afsDet, '.', 'Color', [1 0 0] )
xlabel('X(t_f)')
ylabel('ФД')
title('Зависимость конечного положения траекторий X(t_f) и ФД')
set(gca, 'YScale', 'log')

%% Оценка вероятности времени выхода методом Монте-Карло
% пока что это оценка какая-то другая
prCount = 200;   % количество точек для графика
h = max(max(outs)) / prCount;
levelOut = [0; h * ones(prCount-1,1)];
levelOut = cumsum(levelOut);
% ind = abs(outs(end,:)) > levelOut;
ind = outs(end,:) > levelOut;
probsMCxT = sum(ind,2) / simCount;

tic, fprintf('Вероятности... ')

probsMC = zeros(size(levelOut));
parfor i = 1:size(levelOut,1)
    ind = outs > levelOut(i);
    probsMC(i) = sum(sum(ind)) / (simCount * size(time,2));
end
toc
% figure(5);plot(probsMC)

%% Расчёт ФД и вероятностей методом стрельбы
shootProbs.af = [];
shootProbs.afControl = [];

% левый край
y0 = [x0];

% правый край
% yf = [0.0004533];
% yfs = linspace(0.000001, .00010, 20);
% yfs = linspace(0.00001, .00024, 20);

yfs = linspace(0.00001, h * (prCount *1.25), 20);

timeM = 0.1: 0.1: 1;
tic, fprintf('Оценка вероятности...')
for yf = yfs
    afsShoot = [];
    afsControl = [];
    
    for tM = timeM
        nu = [-120; 120];
        tLocal = 0 : 0.01 : tM;
        
        [shoot] = pathViaShoot(tLocal, sys, nu, y0, yf);
        
        % Оценка вероятности выхода
        mI = size(shoot.syspath(end,:,1),2);   % весь промежуток
        
        af = shoot.control(1:mI,end)' * shoot.control(1:mI,end);  %  * (tLocal(2) - tLocal(1))
        afsControl = [afsControl af];
        
        af = getAF1D(sys, squeeze(shoot.syspath(end,:,:))');  % / (timeM(2) - timeM(1));
        afsShoot = [afsShoot af];
    end
    
    % оценка вероятности времени
    [af, I] = min(afsShoot);
    probExTime = exp(- af / sys.B^2);
    shootProbs.af = [shootProbs.af; probExTime];
 
    [af, I] = min(afsControl);
    probExTime = exp( - af);
    shootProbs.afControl = [shootProbs.afControl; probExTime];
end
toc


%% Рисунки

% Вероятности P(t<T) для заданного уровня
figure(2), clf
subplot(121), hold on, grid on
plot(levelOut, probsMC)
plot(yfs, shootProbs.af, '-', 'Color', [0 .8 .1], 'Marker', 's')
plot(yfs, shootProbs.afControl, '-', 'Color', [.9 0 .1], 'Marker', '*')
plot(levelOut, probsMCxT)
title('Вероятности P(t<T) для заданного X_G')
xlabel('X_G'), ylabel('Частота (МК); Вероятность')

subplot(122), hold on, grid on
plot(levelOut, probsMC)
plot(yfs, shootProbs.af, '-', 'Color', [0 .8 .1], 'Marker', 's')
plot(yfs, shootProbs.afControl, '-', 'Color', [.9 0 .1], 'Marker', '*')
plot(levelOut, probsMCxT)
set(gca, 'YScale', 'log')
xlabel('X_G'), ylabel('Частота (МК); Вероятность')
% legend('Монте-Карло','\intx-Ax','\int u^2', 'Монте-Карло, x(T)', 'Location', 'south')
legend('P^{(1)}_G','\intx-Ax','\int u^2', 'P^{(2)}_G', 'Location', 'southwest')


% Функционалы действия для последнего расчёта методом стрельбы
figure(9), clf, hold on, grid on
plot(timeM, afsShoot)
plot(timeM, afsControl)
xlabel('Время, с'), ylabel('ФД')
set(gca, 'YScale', 'log')
title('Функционалы действия для последнего расчёта методом стрельбы')
legend('\int x-Ax','\int u^2')


% траектория среди всего множества
figure(1)
plot(tLocal, shoot.syspath(end,:))

