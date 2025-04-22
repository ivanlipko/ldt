% ������ ����� � ���������� ��� ��������� ����� � ��������� ���� 
% 
% 2024.02.20 ��������� sys.StateUnit; ����� ���� �������, ����� ���
% �������� � ������ ���������
% 2019.06.26 ������ ������
% 
function [sys, swSpectra, swFreqs] = getKatWindModel(swWindSpeed, yaw, cols)
% swWindSpeed, �/�
% yaw, �������

if nargin == 2
    cols = [6, 5, 3, 1, 4, 2];  % roll first
    % cols = [4, 2, 6, 5, 3, 1];  % pitch first
end
swFreqCount = 500;

[swSpectra, swFreqs, swHs] = waveSpectra(swWindSpeed, swFreqCount);
% disp(['�������� ������ ����� (�) ' num2str(swHs)])

swConstWaveHeight = 5*10^7; % ����������� � �������� ������ ����� Hs
swG2Force = swForceTf(swConstWaveHeight, swSpectra, swFreqs);
swG2Force = ss(swG2Force);
swG2Force.StateName = {'dF','F'};
swG2Force.y = {'Moment'};  % ��� ������������� �����, ��� ��� ������
swG2Force.Notes = '��������� �� �������� ������, � ����� - ����';

% 1000, 100000 - wind speed 10. ��� ��������� ������ �������� ��� ������� �������� �����.
swG2Moment = swForceTf(7*10^4, swSpectra, swFreqs);
swG2Moment.Name = 'Wave to Moment';
swG2Moment.OutputName = 'Moment';
swG2Moment = ss(swG2Moment);
swG2Moment.StateName = {'dM','M'};
swG2Moment.y = {'Force'};  % ��� ������������� �����, ��� ��� ����
swG2Moment.Notes = '��������� �� �������� ���������, � ����� - ������';

% ������ ����������
[~, katFb] = getKatamaran(cols, yaw);

% ���������� ������ ����� � ���������� � ���� �����
sys = connect(swG2Moment, swG2Force, katFb, swG2Moment.InputName, katFb.OutputName);

% sys.C = eye(size(sys.A));

% ������ ������ ���������� ����� (���� � �������� �����) 
cols = [5, 6, 7, 8, 9, 10, 1, 2, 3, 4]; 
sys.A = sys.A(cols, cols);
sys.B = sys.B(cols, :);
sys.C = sys.C(:, cols);
% sys.D = sys.D(cols, :);
sys.StateName = sys.StateName(cols);
sys.StateUnit = sys.StateUnit(cols);

% ������� � ������ ������� (!!! ����������� �� ������ �������������)
sys1 = ss(sys.A, sys.B, eye(size(sys.A,1)), []);
sys1.Name = 'KatWindModel';
sys1.StateName = sys.StateName;
sys1.InputName = sys.InputName;
sys1.OutputName = sys.StateName;
sys1.StateUnit = sys.StateUnit;

sys = sys1;

%%
% A = [ katFb.A    zeros(size(katFb.A,1), 1)  katFb.B * swG2Moment.C'   zeros(size(katFb.A,1), 1)   katFb.B * swG2Force.C'; 
%       zeros(size(swG2Moment.A,1), size(katFb.A,1))  swG2Moment.A   zeros(size(swG2Force.A,1));
%       zeros(size(swG2Moment.A,1), size(katFb.A,1))  zeros(size(swG2Moment.A,1)) swG2Force.A;
%       ];
% B = [swG2Moment.B;
%     swG2Force.B;
%     zeros(size(katFb.A,1), 1) ];
% C = [eye(size(katFb.A))  zeros(size(katFb.A,1), size(swG2Moment.A,1))  zeros(size(katFb.A,1), size(swG2Moment.A,1))];
% sys = ss(A,B,C,[]);
% 
% simOut = lsim(sys, w, time);
% 
% disp(['������� ���������� ' num2str(sum(sum(katSim - simOut))  )])% ������ ���� ������ � ����
% 
% figure(1), clf, hold on
% plot(time, katSim(:,1))
% plot(time, simOut(:,1))