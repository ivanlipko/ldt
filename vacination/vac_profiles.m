% Строим профили критических состояний для модели вакцинации

clear
clc
close all

%% Добавляю пути к скриптам
[path_main, ~, ~] = fileparts(mfilename('fullpath'));
cd(path_main)
addpath(genpath( './../../action_functional_modules' ))

%% Модель 

sys = getVacModel1();

disp(sys.StateName)
