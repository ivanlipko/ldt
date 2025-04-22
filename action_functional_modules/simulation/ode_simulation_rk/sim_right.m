%% Интегрирование с любым временным шагом при наличии внешних возмущений или управления
% Автор: Липко И.Ю. ivanlipko@yandex.ru
%
% 30.03.2021
% 
% 
% Для dt = 0.01; Tf = 2;
% 
% Не используй глобальные переменные, а тем более глобальные классы-объекты
% КОгда использовал глобальную ss-систему время выполнения 0.9-1.1 сек, 
% а при передаче тех же матриц параметрами функций - 0.25-0.27 сек.
% 
% Лучше передавать матрицы, а не объекты. Таким образом уменшил время
% выполнения до около 0.011-0.014 сек.
% 
% Просто добавление глобальных переменных -- увеличило время до 
% 0.015-0.018 сек
% 
% Добавил ещё возмущение по индексу массива -- время выполнения 
% 0.013-0.018 сек
% 
% При частом выполнение наблюдается уменьшение времени вплоть до 0.005 сек.
% Возможно это кеширование.
% 

clear
clc
% close all


%%  Формирование модели объекта разомкнутой системы ДПТ
% Параметры двигателя
dptJ = 0.02;
dptR = 2.0;
dptL = 0.1;
dptkf = 0.075;
dptPsi = 0.7;

disp('Модель объекта разомкнутой системы ДПТ')

dptA = [-dptkf/dptJ      dptPsi/dptJ;
    -dptPsi/dptL     -dptR/dptL];
dptB = [0      -1/dptJ;
    1/dptL    0];
dptC = [1 0];
dptD = [0];

% x' = Ax + Bu; y = Cx + Du.
ssDpt = ss(dptA, dptB, dptC, dptD, 'StateName', {'dc_velocity'; 'dc_I'});
ssDpt.Inputname = {'dc_I'; 'dc_peterb'};
ssDpt.Outputname = 'dc_vel';
ssDpt


%% Синтез LQR оптимального ПИД-регулятора
% Добавление интеграла (Psi * 1/s)
dptExt = [1; tf(1,[1 0])] * ssDpt(1);
dptExt.StateName{1} = 'dc_angle';
% ssDptExtended.Inputname = {'dc_I'; 'dc_peterb'};
dptExt.Outputname = {'dc_vel'; 'dc_angle'};
dptExt

% если что - смотреть лаб №4. 10 семестр
q = 20.0;   % параметры LQR критерия
ro = 0.01;
Q = [1 0;
    0 q];
R = [ro];
[K, S, e] = lqry(dptExt, Q, R);

dptCPid = ss( K * [ tf(1,[1 0]) 0 0;  0 1 0;  0 0 1 ] );
% Дополнение выхода разомкнутой системы ДПТ вектором состояния
dptAug = augstate( ssDpt );
dptAug = dptAug * append( dptCPid, 1 );

% Формирование закнутой системы
dptAugFb = feedback( dptAug, eye(3), [1 2 3], 1:3 );
dptAugFb = dptAugFb( 1, [1 4] );

% для x' = Ax + Bu + Gw; (w - noise)
% значения взял от балды
dptAugFbG = [
    0.1 0;
    0 0.1
    0 0];


%% Моделирование
dt = 0.01;
Tf = 2;
time = 0 : dt : Tf;
% ptb = dptAugFb.B* wgn(2, length(time), 1);
ptb = dptAugFbG * zeros(2, length(time));
% ptb = zeros(2, length(time));   % если понадобится изменять матрицу G в момент моделирваония

cntrl = .1*dptAugFb.B* wgn(2, length(time), 1);  
% cntrl = wgn(2, length(time), 1);   % если понадобится изменять матрицу G в момент моделирваония

x0 = [1; 5; 3];    

state = x0;
% out = [state'; zeros(length(time), length(state)) ];
out = zeros(length(time), length(state));

right = @rightPeterbFullSystem;

sysA = dptAugFb.A;
sysB = dptAugFb.B;
sysG = dptAugFbG;

pbc = length(time) / 10;
pbs = 0; pb = 0;

i = 1;

% profile on
tic
for t = time
    pbs = pbs + 1;
    if pbs > pbc
        fprintf('%d.', pb);
        pbs = 0;
        pb = pb + 1;
    end
    
    out(i,:) = state;    
    state = rkStepPeterbFullSystem(t, dt, state, right, sysA, sysB, sysG, cntrl(:,i), ptb(:,i));

    i = i + 1;
end
% time(end+1) = Tf+dt;
% out(i,:) = state;

toc

% profile viewer

figure(1), clf, grid on, hold on
plot(time, out)

%% проверка с lsim: сошлось
% time(end) = [];
% outlsim = lsim(dptAugFb, cntrl, time, x0);
% outlsim = lsim(dptAugFb, zeros(length(time),2), time, x0);
% plot(time, outlsim)


%%
[path2module, ~, ~] = fileparts(mfilename('fullpath'));
cd(path2module)
addpath(genpath( '../../action_functional_modules' ))
clear path2module
% restoredefaultpath 

rng(124)
simCount = 5;

ptb = 20*randn(length(time), size(dptAugFb.B,2), simCount);
outs = simulate_model_parfor_right(dptAugFb, dptAugFbG, x0, time, ptb, @rightdpt, simCount);

figure(2), clf, grid on, hold on
plot(time, squeeze(outs(:,2,:)))


%%
function dx = rightdpt(t, x, A,ptb)
% dx = A * x + B*cntrl + G*ptb;

dx = A * x + ptb;
end