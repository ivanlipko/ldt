%% ������ ������-������ ����� � ���������
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
% �� ����
% x=[x_1  x_2 ]^T � ������ ���������, ���������� ����� � �������� ��������������, 
% ? � ����������� ����������� �����, ? � ����������� ����� ��������, 
% ? � ����������� �������� �����, ? � ����������� ��������������� ����� ��������, 
% g_1  � g_2 � ��������� ���������� �����, ?(t) � ����� ��� (���������� �����).

function [right, ps, right_migrt] = getPreyPredator_nonlin(al, bt, gm, dl, dtl)
% ������ �����
ps = [gm/dl al/bt];

% ���������� ������
right = @(t,x)[al*x(1) - bt*x(1)*x(2);    % ���������� �����
           -gm*x(2) + dl*x(1)*x(2)] ;  % ���������� ��������

% from art. Numerical comp. p 11
% right = @(t,x)[al*x(1) - bt*x(1)*x(2) + dtl;    % ���������� �����
%               -gm*x(2) + dl*x(1)*x(2) + dtl] ;  % ���������� ��������

       
right_migrt = @right_migration;

end
%%

% ���������� ������ � �������� ��� ����
% function dx = right_migration(t,x,w,wtime, al,bt,gm,dl,g1,g2)
% f = interp1(wtime,w,t, 'previous'); % Interpolate the data set (ft,f) at time t
% dx = zeros(2,1);
% dx(1) =  al*x(1) - bt*x(1)*x(2) + g1*f(1);    % ���������� �����
% dx(2) = -gm*x(2) + dl*x(1)*x(2) + g2*f(2);         % ���������� ��������
% end

% function dx = right_migration(t,x,p)
% al = p.al;
% bt = p.bt;
% gm = p.gm;
% dl = p.dl;
% g1 = p.g1;
% g2 = p.g2;
% w = p.w;

% function dx = right_migration(t,x,  al,bt,gm,dl,g1,g2, w)
% dx = zeros(2,1);
% dx(1) =  al*x(1) - bt*x(1)*x(2) + g1*w(1);    % ���������� �����
% dx(2) = -gm*x(2) + dl*x(1)*x(2) + g2*w(2);         % ���������� ��������
% end

% from art. Numerical comp. p 11
function dx = right_migration(t,x,  al,bt,gm,dl,g1,g2, w, dtl)
dx = zeros(2,1);
dx(1) =  al*x(1) - bt*x(1)*x(2) + dtl + g1*w(1);    % ���������� �����
dx(2) = -gm*x(2) + dl*x(1)*x(2) + dtl + g2*w(2);         % ���������� ��������
end
