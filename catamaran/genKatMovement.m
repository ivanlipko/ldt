% Генерация движений катамарана
% 
% Автор:
%   Липко Иван Юрьевич
%   ivanlipko@yandex.ru
% 
% TODO:
% 
% История версий:
%   2021.07.05 ДОбавлено: выходные параметры ptrbSim, w
%              Почистил немного код
%   10.10.2018 Разделение на части
%   03.10.2018 Начальная версия.
% 
% function [katFb, katSim] = genKatMovement(time, swWindSpeed)  2018
% function [katFb, katSim, w] = genKatMovement(time, swWindSpeed)

function [katFb, katSim, ptrbSim, w] = genKatMovement(time, swWindSpeed)

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
[~, katFb] = getKatamaran(cols, 0);

%% Симулирование
w = wgn(size(time,2), 1, 0);
% Внешние возмущеня. Берём из файла чтобы была повторяемость
% load('./data/sim_3-10-2018 11.23.15.mat', 'w');  % time is 50*10^4
% load('./data/sim_3-10-2018 12.40.54.mat', 'w'); % time is 1*10^3
% swForce и swMoment не повторяются для одного и того же шума!

swForce = lsim(swSpec2Force, w, time);
swMoment = lsim(swSpec2Moment, w, time);

% реакция объекта с управлением
% ptrb = [swForce swMoment]; % возмущение для cols == [1,2,3,4]
ptrbSim = [swMoment swForce]; % возмущение для cols == [4,2,1,3]
katSim = lsim(katFb, ptrbSim, time);

% реакция объекта без управления
% disturb = [zeros(length(time),1) zeros(length(time),1) swForce swMoment];
% disturb = disturb(:, cols); % входной сигнал надо переставлять
% katmaranSimOut_nofb = lsim(katamaran, disturb, time);

% clear swForce swMoment



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