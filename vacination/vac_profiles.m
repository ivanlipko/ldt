% ������ ������� ����������� ��������� ��� ������ ����������

clear
clc
close all

%% �������� ���� � ��������
[path_main, ~, ~] = fileparts(mfilename('fullpath'));
cd(path_main)
addpath(genpath( './../../action_functional_modules' ))

%% ������ 

sys = getVacModel1();

disp(sys.StateName)
