% ������������� ������ ������� �����-�����
%
% ������ ����������� ������� ������� ������ �� �������
% P (t < T), t = {t: x > XThr}
% 
%  (!!!) ���� ���� ����������, ����� ������� ������� �������. ��� ������
%  ��������, ����� ���� ������������, ������������ ��������� � ������
%  �����-����� � � ������ ��������

[path2module, ~, ~] = fileparts(mfilename('fullpath'));
cd(path2module)
addpath(genpath( '../../action_functional_modules' ))

%% �������������. ������
time = 0 : 0.01 : 1;
simCount = 5000;

x0 = [0];
sys = get1Dsys();

% [outs, w] = simulate_model(sys, x0, simCount, time, 3);
[outs, w] = simulate_model_parfor(sys, time, x0, simCount);

% outDet = lsim(sys, zeros(size(time,2),1), time, x0);   % ������������� ��������

outs = squeeze(outs);
w = squeeze(w);


%% �� S = int(dot_x - b(x)) � �����������. ������
tic, fprintf('��... ')

afsMC =  zeros(simCount,1);
parfor i = 1:simCount
    afsMC(i) = getAF1D(sys, outs(:,i));
end
toc

% probs = exp(- afsMC / sys.B^2);


%% ����������. �������
tic, fprintf('Plot...')

figure(1), clf, hold on, grid on
plot(time, outs, 'Color', [.5 .5 .5 .1])
title('����������'), xlabel('�����, �')

toc

%% ����������� ��������� ������ ���������� � �� ������� �����-�����
figure(4), clf, hold on, grid on
plot(outs(end,:), afsMC, '.', 'Color', [.5 .5 .5 .5])
% plot(outDet(end), afsDet, '.', 'Color', [1 0 0] )
xlabel('X(t_f)')
ylabel('��')
title('����������� ��������� ��������� ���������� X(t_f) � ��')
set(gca, 'YScale', 'log')

%% ������ ����������� ������� ������ ������� �����-�����
% ���� ��� ��� ������ �����-�� ������
prCount = 200;   % ���������� ����� ��� �������
h = max(max(outs)) / prCount;
levelOut = [0; h * ones(prCount-1,1)];
levelOut = cumsum(levelOut);
% ind = abs(outs(end,:)) > levelOut;
ind = outs(end,:) > levelOut;
probsMCxT = sum(ind,2) / simCount;

tic, fprintf('�����������... ')

probsMC = zeros(size(levelOut));
parfor i = 1:size(levelOut,1)
    ind = outs > levelOut(i);
    probsMC(i) = sum(sum(ind)) / (simCount * size(time,2));
end
toc
% figure(5);plot(probsMC)

%% ������ �� � ������������ ������� ��������
shootProbs.af = [];
shootProbs.afControl = [];

% ����� ����
y0 = [x0];

% ������ ����
% yf = [0.0004533];
% yfs = linspace(0.000001, .00010, 20);
% yfs = linspace(0.00001, .00024, 20);

yfs = linspace(0.00001, h * (prCount *1.25), 20);

timeM = 0.1: 0.1: 1;
tic, fprintf('������ �����������...')
for yf = yfs
    afsShoot = [];
    afsControl = [];
    
    for tM = timeM
        nu = [-120; 120];
        tLocal = 0 : 0.01 : tM;
        
        [shoot] = pathViaShoot(tLocal, sys, nu, y0, yf);
        
        % ������ ����������� ������
        mI = size(shoot.syspath(end,:,1),2);   % ���� ����������
        
        af = shoot.control(1:mI,end)' * shoot.control(1:mI,end);  %  * (tLocal(2) - tLocal(1))
        afsControl = [afsControl af];
        
        af = getAF1D(sys, squeeze(shoot.syspath(end,:,:))');  % / (timeM(2) - timeM(1));
        afsShoot = [afsShoot af];
    end
    
    % ������ ����������� �������
    [af, I] = min(afsShoot);
    probExTime = exp(- af / sys.B^2);
    shootProbs.af = [shootProbs.af; probExTime];
 
    [af, I] = min(afsControl);
    probExTime = exp( - af);
    shootProbs.afControl = [shootProbs.afControl; probExTime];
end
toc


%% �������

% ����������� P(t<T) ��� ��������� ������
figure(2), clf
subplot(121), hold on, grid on
plot(levelOut, probsMC)
plot(yfs, shootProbs.af, '-', 'Color', [0 .8 .1], 'Marker', 's')
plot(yfs, shootProbs.afControl, '-', 'Color', [.9 0 .1], 'Marker', '*')
plot(levelOut, probsMCxT)
title('����������� P(t<T) ��� ��������� X_G')
xlabel('X_G'), ylabel('������� (��); �����������')

subplot(122), hold on, grid on
plot(levelOut, probsMC)
plot(yfs, shootProbs.af, '-', 'Color', [0 .8 .1], 'Marker', 's')
plot(yfs, shootProbs.afControl, '-', 'Color', [.9 0 .1], 'Marker', '*')
plot(levelOut, probsMCxT)
set(gca, 'YScale', 'log')
xlabel('X_G'), ylabel('������� (��); �����������')
% legend('�����-�����','\intx-Ax','\int u^2', '�����-�����, x(T)', 'Location', 'south')
legend('P^{(1)}_G','\intx-Ax','\int u^2', 'P^{(2)}_G', 'Location', 'southwest')


% ����������� �������� ��� ���������� ������� ������� ��������
figure(9), clf, hold on, grid on
plot(timeM, afsShoot)
plot(timeM, afsControl)
xlabel('�����, �'), ylabel('��')
set(gca, 'YScale', 'log')
title('����������� �������� ��� ���������� ������� ������� ��������')
legend('\int x-Ax','\int u^2')


% ���������� ����� ����� ���������
figure(1)
plot(tLocal, shoot.syspath(end,:))

