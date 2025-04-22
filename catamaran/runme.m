% ������ ������� ����������� ���������

% clear
% clc
% close all

%% �������� ���� � ��������
[path_main, ~, ~] = fileparts(mfilename('fullpath'));
cd(path_main)
addpath(genpath( './../../action_functional_modules' ))


%% ������ ����� � ���������� � ����� �������
% 
%  ����� �������� ����� ������������ ��������� ����� �������� �������
%  ������� ���������� ������ �������
% 
arWindSpeed = 10;
sys = getKatWindModel(arWindSpeed, 0);

disp(sys.StateName)

if 1
    return
end
%% ��������� ���������� ������ (���������� + ���������)
time = 0:0.1:50;
[~, katSim, ~, w] = genKatMovement(time, arWindSpeed);
simOut = lsim(sys, w, time);

disp(['������� ���������� ' num2str(sum(sum(katSim - simOut))  )])  % ������ ���� ������ � ����

% figure(1), clf, hold on
% plot(time, katSim(:,1))
% plot(time, simOut(:,1))

% clearvars time simOut katOut w


%% ������� �������� (� ������������ �������� ��� ���������� ������)
xf = 5;  % 3:1:15;  % ��������� ��������
% swWindSpeed = 1:3:10;
arWindSpeed = 1:3:25;

Tf = 10;  % seconds
dT = 0.1;

f = figure(1); f.Position = [300 350 507 527]; clf
subplot(211), hold on, grid on
subplot(212), hold on, grid on

f = figure(2); f.Position = [450 350 507 527]; clf
subplot(211), hold on, grid on
subplot(212), hold on, grid on

% f = figure(3); %f.Position = [450 350 507 527]; clf
% clf, hold on

% for xTf = xTf
for wS = arWindSpeed
    [sys, swSpectra, swFreqs] = getKatWindModel(wS, 45);
%     [pfl.time, pfl.path, pfl.peterb, pfl.af] = getProfile(sys, sys.B, xf, Tf, dT);
    pfl = getProfile(sys, sys.B, xf, Tf, dT);
    
    figure(1)
    subplot(211)
    plot(pfl.time, pfl.apath(1,:))
    subplot(212)
    plot(pfl.time, pfl.apath(5,:))
    
    figure(2)
    subplot(211)
    plot(pfl.time, pfl.apath(2,:))
    subplot(212)
    plot(pfl.time, pfl.apath(6,:))
    
%     subplot(211)
%     plot(pfl.apath(1,:), pfl.apath(2,:))
%     subplot(212)
%     plot(pfl.apath(5,:), pfl.apath(6,:))
    
%     figure(3)
%     plot(swFreqs, swSpectra, 'Linewidth', 2)
end

figure(1)
subplot(211),
xlabel('�����, �'), ylabel(sys.StateName(1))
% xlabel('�����, �'), ylabel('����, ^o')
% xlabel('�����, �'), ylabel('���� ����������, ^o')
legend( string(num2cell(arWindSpeed)) + ' �/�', 'Location', 'West')
subplot(212)
xlabel('�����, �'), ylabel(sys.StateName(5))
% xlabel('�����, �'), ylabel('���� ����������, ^o')
% xlabel('�����, �'), ylabel('����, ^o')
legend( string(num2cell(arWindSpeed)) + ' �/�', 'Location', 'West')

figure(2)
subplot(211),
xlabel('�����, �'), ylabel(sys.StateName(2))
% xlabel('�����, �'), ylabel('�������� �����, ^o/�')
% xlabel('�����, �'), ylabel('�������� ����������, ^o/�')
legend( string(num2cell(arWindSpeed)) + ' �/�', 'Location', 'West')
subplot(212)
xlabel('�����, �'), ylabel(sys.StateName(6))
% xlabel('�����, �'), ylabel('�������� ����������, ^o/�')
% xlabel('�����, �'), ylabel('�������� �����, ^o/�')
legend( string(num2cell(arWindSpeed)) + ' �/�', 'Location', 'West')

% figure(3)
% title('������ ��������'), grid on
% xlabel('�������, ���/�'), ylabel('S(\omega), ���/�')
% legend( string(num2cell(swWindSpeed)) + ' �/�', 'Location', 'East')
% axis([0 1.5 0 50])
% axis auto



%% ������� �������� (��������� ����� ��� ���������� ������)
arxf = 5:2:10;  % ��������� ��������
arWindSpeed = 1:3:10;
% swWindSpeed = 10:5:25;

Tf = 10;  % seconds
dT = 0.1;

f1 = figure(1); f1.Position = [300 350 600 350]; clf
subplot(211), hold on, grid on, g1 = gca;
subplot(212), hold on, grid on, g2 = gca;

f = figure(2); f.Position = [450 350 600 350]; clf
subplot(211), hold on, grid on, g3 = gca;
subplot(212), hold on, grid on, g4 = gca;

f = figure(3); f.Position = [450 350 600 350]; clf
hold on, grid on, g5 = gca;

for xf = arxf
    for wS = arWindSpeed
        [sys, swSpectra, swFreqs] = getKatWindModel(wS, 90);
%         [pfl.time, pfl.path, pfl.peterb, pfl.af] = getProfile(sys, sys.B, xf, Tf, dT);
        pfl = getProfile(sys, sys.B, xf, Tf, dT);

        figure(1)
        subplot(211)
        plot(pfl.time, pfl.apath(1,:))
        subplot(212)
        plot(pfl.time, pfl.apath(5,:))

        figure(2)
        subplot(211)
        plot(pfl.time, pfl.apath(2,:))
        subplot(212)
        plot(pfl.time, pfl.apath(6,:))
        
        figure(3)
        plot(pfl.time, normalize(pfl.apath(1,:), 1, 'norm'), 'Color', [.3 .3 .3] )
        plot(pfl.time, normalize(pfl.apath(8,:), 1, 'norm'), 'Color', [.0 .8 .3] )

    end
    g1.ColorOrderIndex = 1;
    g2.ColorOrderIndex = 1;
    g3.ColorOrderIndex = 1;
    g4.ColorOrderIndex = 1;
    g5.ColorOrderIndex = 1;
end

figure(1)
subplot(211),
xlabel('�����, �'), ylabel('����, ^o')
% xlabel('�����, �'), ylabel('���� ����������, ^o')
legend( string(num2cell(arWindSpeed)) + ' �/�', 'Location', 'SouthWest')

subplot(212)
xlabel('�����, �'), ylabel('���� ����������, ^o')
% xlabel('�����, �'), ylabel('����, ^o')
legend( string(num2cell(arWindSpeed)) + ' �/�', 'Location', 'SouthWest')

figure(2)
subplot(211),
xlabel('�����, �'), ylabel('�������� �����, ^o/�')
% xlabel('�����, �'), ylabel('�������� ����������, ^o/�')
legend( string(num2cell(arWindSpeed)) + ' �/�', 'Location', 'NorthWest')
subplot(212)
xlabel('�����, �'), ylabel('�������� ����������, ^o/�')
% xlabel('�����, �'), ylabel('�������� �����, ^o/�')
legend( string(num2cell(arWindSpeed)) + ' �/�', 'Location', 'SouthWest')

figure(3)
xlabel('�����, �')
legend( '���� ����������', '����', 'Location', 'NorthWest')


%%

%% ������� ���������� � ��������
arxf = 5;  % ��������� ��������
arWindSpeed = 1:5:10;
% swWindSpeed = 10:5:20;
% swWindSpeed = [1  25];

Tf = 40;  % seconds
dT = 0.1;

figure(4), clf, hold on, grid on, g7 = gca;

for xf = arxf
    for wS = arWindSpeed
        [sys, swSpectra, swFreqs] = getKatWindModel(wS, 0);
%         [pfl.time, pfl.path, pfl.peterb, pfl.af] = getProfile(sys, sys.B, xf, Tf, dT);
        pfl = getProfile(sys, sys.B, xf, Tf, dT);
        figure(4)
        plot(pfl.time, normalize(pfl.apath(1,:), 1, 'norm'), 'Color', [.3 .3 .3] )
        plot(pfl.time, normalize(pfl.apath(8,:), 1, 'norm'), 'Color', [0 .8 .3]  )
        plot(pfl.time, normalize(pfl.ptb, 1, 'norm'), 'Color', [.8 .3 .3]  )

%   �������� � �������� �� ����� - �� ���� ������� �������� - ���� ���������, ���-�� ����� �� ���      
%         
%         pfl = getProfileFrom(3, sys.A, sys.B, [xf; zeros(size(sys.A,1)-1, 1)], Tf, dT);
%         figure(4)
% %         plot(pfl.time, normalize(pfl.path(1,:), 1)', 'Color', [.3 .3 .3] )
% %         plot(pfl.time, normalize(pfl.path(8,:), 1)', 'Color', [0 .8 .3]  )
% %         plot(pfl.time, normalize(pfl.ptb, 1, 'norm'), 'Color', [.8 .3 .3]  )
%         plot(pfl.time, pfl.path(1,:), 'Color', [.3 .3 .3] )
%         plot(pfl.time, pfl.path(8,:), 'Color', [0 .8 .3]  )


    end
%     g1.ColorOrderIndex = 1;

end

xlabel('�����, �')
legend( '���� ����������', '����', '���������� ����', 'Location', 'NorthWest')
pfl.ptb * pfl.ptb'
pfl.apath(8,:) * pfl.apath(8,:)'



%%

%% �-������� �������� ������ �� "��������" �-������� �������� ��������
%  � ��� ��� ������������ ��� �������� �����������? ���� ���������� ����
%  ����� �-�������, � ����-�-���� �� ���.
pfl = getProfile(sys, sys.B, 9, Tf, dT/10);
pfl1 = getProfile(sys, sys.B, 7, Tf, dT/10);
pfl1 = getProfile(sys, sys.B, 7, Tf, dT/10);

figure(5), clf, grid on, hold on, title('�-������� ������ ������� �� ������� �� �������')
plot(pfl.time, pfl.apath(1,:))
plot(pfl1.time, pfl1.apath(1,:))
plot(pfl1.time-.85, pfl1.apath(1,:))

figure(6), clf, grid on, hold on, title('�-������� � ��������� ���������. �����������')
plot(pfl.apath(1,:), pfl.apath(2,:))
plot(pfl1.apath(1,:), pfl1.apath(2,:))
plot(pfl1.apath(1,:)*1.3, pfl1.apath(2,:)*1.3)

figure(7), clf, grid on, hold on, title('���, ��� �������� � �����������?')
plot(pfl.time, pfl.aaf)
plot(pfl1.time, pfl1.aaf)

% out = pfl.aaf - pfl.aaf(3972);  ������� ��������� �� �� ���������� � �����
out = pfl.aaf;
out = [out; pfl1.aaf];
out = log(log(out+1));
figure(8), clf, grid on, hold on, title('���, ��� �������� � �����������?')
scatter(pfl.apath(1,:), pfl.apath(2,:), [], out(1,:), 'filled')
scatter(pfl1.apath(1,:), pfl1.apath(2,:), [], out(2,:), 'filled')
line([7 7],[0 15])
colorbar



%% -----------------------------------------------------------------------
%% �-������� � ������������� �����-�����.
% 
%% -----------------------------------------------------------------------
%  ������������� ���

x0 = zeros(1, size(sys.A, 1));
thr = 4;
% statevec = 1 : size(sys.A, 1);
statevec = 1 : 6;

rng(123)

mmcCount = 300; batchSize = 500;
% exceedCount = 300;  batchSize = 500;  tFinal = 600;
% ����� ������������� 3000 Elapsed time is 101.242013 seconds.

tFinal = 600;
dT = 0.01;
simTime = 0 : dT : tFinal;

% ��������� �� ������� ������������
tic, fprintf('Simulation... ')
[mmc_outs, ~, mmc_w] = simulate_model_parforNd_thr(sys, simTime, x0, thr, statevec, mmcCount, batchSize, @simfun);
mmcCount = size(mmc_outs,3);  % update count
toc

% ����������, ������� ����������
figure(21), clf, hold on, grid on, title('����������, ������� ���������� � ��� (������)')
xlabel(sys.StateName(1)), ylabel(sys.StateName(2))
plot( squeeze( mmc_outs(:,1,1:10)) )

% %% ����� ����������� ������������� � ��� ����������� ������� ��������� 
% % �� 10, � ����� 6 - ��� ����� ������� ���, ������� ���� ���������������� ������
% for i=1:mmcCount
%     simOut = lsim(sys, mmc_w(:,1,i), simTime);
% end

%% �������� ������� � �����������
exdSizeTime = 15;
indent = 1500;
exdTime = 0:dT:exdSizeTime-dT;

exdPathes = [];
exdPtb = [];

for i=1:mmcCount
%     [exd, timesAbove, isAbove] = findAbove_v2(sys, mmc_outs(:,:,i), simTime, winSizeTime, thr, indent);
    [exd, timesAbove, isAbove] = findAbove_v2(size(sys.C,1)+1, cat(2, mmc_outs(:,:,i), mmc_w(:,:,i)), simTime, exdSizeTime, thr, indent);
    exdPathes = cat(3, exdPathes, exd);
end

% clear winSizeInds 

%%
% figure(22), clf, hold on, grid on
% xlabel(sys.StateName(1)), ylabel(sys.StateName(2))
% plot( squeeze( exceed_outs(:,1,1:10)) )


% %% ������ �� ����� ��������� ���������� �� ������ � �� � ��������� � �-��������
pfl.time = 0:dT:15;
pfl = getProfile_v2(sys, pfl.time, thr);

figure(22), clf
subplot(211), hold on, grid on, title(['��������� � ��������� � ���, ����� ', num2str(arWindSpeed)])
pl1 = line([0 exdSizeTime],[thr thr], 'Linestyle', '--', 'Color',[1 0 0]);
pl2 = plot(exdTime, exdPathes(:,1,1), 'Color',[.5 .5 .5 .1]);

for i=2:size(exdPathes,3)
    plot(exdTime, exdPathes(:,1,i), 'Color',[.5 .5 .5 .1])
end

pl3 = plot(pfl.time, pfl.apath(1,:), 'Color',[.2 .3 .8], 'Linewidth', 2);
pl4 = plot(exdTime, mean(exdPathes(:,1,:),3), 'Color',[.3 .3 .4], 'Linestyle', '-', 'Linewidth', 2);

xlabel('�����', 'FontSize', 16), ylabel(sys.StateName(1), 'FontSize', 16)
legend([pl1, pl2, pl3, pl4], '������� ��', '���', '���������', '������� ���', 'Location', 'South')

% ������ ������ ���� - ���������� �� ������ � �� � ��������� � �-��������
subplot(212), hold on, grid on, title(['����� ��� � ��������� � ���, ����� ', num2str(arWindSpeed)])
% pl1 = line([0 winSizeTime],[thr thr], 'Linestyle', '--', 'Color',[1 0 0]);
pl2 = plot(exdTime, exdPathes(:,7,1), 'Color',[.5 .5 .5 .1]);

for i=2:size(exdPathes,3)
    plot(exdTime, exdPathes(:,7,i), 'Color',[.5 .5 .5 .1])
end

pl3 = plot(pfl.time, pfl.ptb, 'Color',[.2 .1 .7]);

xlabel('�����', 'FontSize', 16), ylabel('\xi', 'FontSize', 16)
legend([pl2, pl3], '���', '���������', 'Location', 'South')


% %% ������ �� ����� ��������� ���������� �� ������ � �� � ��������� � �-�������� (���������)
figure(23), clf, hold on, grid on, title(['��������� � ��������� � ���, �����', num2str(arWindSpeed)])
% pl1 = line3([0 winSizeTime],[thr thr], 'Linestyle', '--', 'Color',[1 0 0]);
pl2 = plot3(exdPathes(:,1,1), exdPathes(:,3,1), exdPathes(:,5,1), 'Color',[.5 .5 .5 .1]);

for i=2:size(exdPathes,3)
    plot3(exdPathes(:,1,i), exdPathes(:,3,i), exdPathes(:,5,i), 'Color',[.5 .5 .5 .1]);
end

pl3 = plot3(pfl.apath(1,:), pfl.apath(3,:), pfl.apath(5,:), 'Color',[.2 .1 .7], 'LineWidth', 2);

xlabel(sys.StateName(1), 'FontSize', 16)
ylabel(sys.StateName(3), 'FontSize', 16)
zlabel(sys.StateName(5), 'FontSize', 16)
legend([pl2, pl3,], '���', '���������')


%% �������� ��-�� ��� ��� � ����������
% af_pfl = calcAFofPath(pfl.apath', sys.A, sys.B, pfl.time);
% calcAFofPtb(pfl.ptb, pfl.time);
% calcAFofPtb(exdPathes(:,7,1)', 0:dT:winSizeTime-dT)

padding = 1000;
padding * dT

af_exds = zeros(1,size(exdPathes,3));

for i=1:size(exdPathes,3)
    af_exds(i) = calcAFofPtb(exdPathes(end-padding:end,7,i)', exdTime(end-padding:end));
    
% af_exds(i) = calcAFofPtb(exdPathes(:,7,i)', 0:dT:winSizeTime-dT);

%     af_exds(i) = calcAFofPath(exdPathes(:,:,i), ...
%         sys.C * sys.A * sys.C', ...
%         sys.C * sys.B, ...
%         pfl.time);
end

figure(24), clf, hold on, title('�������������� ���������� � ��������� � ���')
histogram(af_exds)
line([pfl.aaf(end-padding) pfl.aaf(end-padding)],[0 50], 'Color', [1 0 0 ])
legend('���','���������')

%
%
%% ������� ������������� ��� ���
function [out, w] = simfun(sys, time, x0)
% w = randn(length(time), size(sys.B,2));
% out = lsim(sys, w, time, x0);
arWindSpeed = 10;
% [~, out, w] = genKatMovement(time, arWindSpeed);
[~, out, ~, w] = genKatMovement(time, arWindSpeed);
end