% Генерация движений катамарана
% 
% Автор:
%   Липко Иван Юрьевич
%   ivanlipko@yandex.ru
% 
% TODO:
%  модуль называется "генерация траектории для катамарана", а получается
%  что я и сам катамаран внутри генерирую. не совсем хорошо наверно
%  получается. надо чтобы катамаран создавался отдельно. Т.е. передавать
%  его в качестве параметра
% 
% История версий:
%   10.10.2018 Разделение на части
%   03.10.2018 Начальная версия.
% 
% dT = .1;
% tFin = 5*10^2;
% swWindSpeed = 10;
% 
% time = 0:dT:tFin;
% 
% 
function [sys, sim, w] = genKatMovement(time, swWindSpeed, angleDegPeterb)
if nargin < 3
    angleDegPeterb = 0;
    warning('angleDegPeterb = 0;')
end

%% Добавляю пути к скриптам
[path_main, ~, ~] = fileparts(mfilename('fullpath'));
cd(path_main)
addpath(genpath( fullfile(path_main, 'data')))
% addpath(genpath( fullfile(path_main, 'rk4Step')))
addpath(genpath( fullfile(path_main, '../shareAlgorithms')))
addpath(genpath( fullfile(path_main, 'disturbances')))

% return

%% Спектр волнения
% swWindSpeed = 10;
swFreqCount = 500;
[swSpectra, swFreqs, ~] = waveSpectra(swWindSpeed, swFreqCount);
sw_const_waveHeight = 5*10^7; % приближение к значащей высоте волны Hs
swSpec2Force = swForceTf(sw_const_waveHeight, swSpectra, swFreqs);
% 1000, 100000 - wind speed 10. Эти параметры должны меняться для каждого
% значения ветра.
swSpec2Moment = swForceTf(7*10^4, swSpectra, swFreqs);

%% Модели катамарана
% global katamaranFb
cols = [6, 5, 3, 1, 4, 2]; % roll first
% sys_state_num = 1;
[~, sys] = getKatamaran(cols, angleDegPeterb);

%% Симулирование
tic

% Внешние возмущеня. Берём из файла чтобы была повторяемость
w = wgn(size(time,2), 1, 0);
% load('./data/sim_3-10-2018 11.23.15.mat', 'w');  % time is 50*10^4
% load('./data/sim_3-10-2018 12.40.54.mat', 'w'); % time is 1*10^3
% swForce и swMoment не повторяются для одного и того же шума!
swForce = lsim(swSpec2Force, w, time);
swMoment = lsim(swSpec2Moment, w, time);

% реакция объекта с управлением
% ptrb = [swForce swMoment]; % возмущение для cols == [1,2,3,4]
ptrb = [swMoment swForce]; % возмущение для cols == [4,2,1,3]
sim = lsim(sys, ptrb, time);

% реакция объекта без управления
% disturb = [zeros(length(time),1) zeros(length(time),1) swForce swMoment];
% disturb = disturb(:, cols); % входной сигнал надо переставлять
% katmaranSimOut_nofb = lsim(katamaran, disturb, time);
toc
clear swForce swMoment

%% Сохранение данных симулирования
% save(['./data/sim_',string(datetime('now', 'Format','d-M-y HH.mm.ss')).char,'.mat'], ...
%   'w', 'katmaranSimOut', 'time')

%% Рисунки:
%% Действующее возмущение
% figure(43), clf
% plot(time, swForce)
% grid on, title('Wave to Force, N'), xlabel('Time, seconds'), ylabel('Force')

%% какой белый шум на входе
% figure(46), clf, hold on
% plot(time, w, 'Color', [1 0 0 ])
% % plot(time, swForce)
end