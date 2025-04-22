%% ������ ������-������ ����� � ���������
% �������� ������� (https://ru.wikipedia.org/wiki/������ ����� � ���������)
% x0 = [ps(1) ps(2)] - x0;
% 
% al = 4;
% bt = .05;
% gm = 2;
% dl = .05;
% 
% al = 1;
% bt = .14;
% gm = .42;
% dl = .07;
%
%% ������� � word:
% https://www.youtube.com/watch?v=3jfLmmcG8fw     ��  10 ������
%             \matrix(0 & -bt*gm/dl @ dl*al/bt & 0)
% ������ ������ �����^               � ����� ����� ^
% 

%%
function [sys, ps] = getPreyPredator_linear(al, bt, gm, dl, g1, g2)
% ������ �����
ps = [gm/dl al/bt];

A = [0 -bt*gm/dl;
     dl*al/bt 0];

G = [ g1  0     
      0   g2];  
C = eye(size(A));
D = zeros(size(G));

sys = ss(A,G,C,D);
sys.Name = '������-������ (�����-���������). ������ � ����������� �� ������������ �����';
sys.StateName = {'������', '�������' };
sys.StateUnit = {'��','��'};
sys.InputName = {'������-��������', '�������-��������'};
sys.OutputName = sys.StateName;
end

