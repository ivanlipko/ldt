% ������ ���������������� ��������� ������� � ������� ���������
%
%  curState - ������ ��������� ������� � ������� ������ ������� (!!!)
%  mu - ������ �������� �����
%     size(mu) == size(sysState)
%
%  ����������:
%   2022.05.31
%       ������������ indMinDist � minDistInd
%   2022.04.28
%       �������� ����������� ������� � ���������� path
%   2020.08.28
%       ������� � �������� �������� ���������, ���������� ���� ���������.
%       �����������
%   2018.10.27
%       �����������. �������������� ������ �� calcH � getDistToProfl
%   2018.04.18
%       ���������. ����������� ������ ������� � ������� mu: �� ���������� �
%       ����������� �� �� ��������. ones(size(kedTime,2), size(mu,2))
%   2018.03.20
%       ��������� ������ ������� � 31 ��� �� ���� ������������� ���������
%       ��������. ������� ����.

function [minDist, minDistTime, minDistInd, v] = getDistToProfl(time, path, cstate, mu)
assert( sum(size(mu) == size(path(1,:)) ) == 2, 'sum(size(mu) == size(prfl(1,:))) == 2 ' )
assert( size(cstate,2) == size(path,2) , 'must be TRUE size(cstate,2) == size(path,2)' )
assert( size(time,2) == size(path,1) , 'must be TRUE size(time,2) == size(path,1) ' )

dist = mu .* ones(size(time,2), size(mu,2)) .* ( (cstate - path).^2 );   % euclid
% dist = mu .* ones(size(time,2), size(mu,2)) .* sum( abs(cstate - path),2 );     % manhatten
% sum(abs(x-y))

% sqrt??? ����� ����� ��� ����� �������
[minDist, minDistInd ] = min(sum(dist, 2));
minDistTime = time(minDistInd);

v.minDist = minDist;
v.minDistTime = minDistTime;
v.indMinDist = minDistInd;
end



%{
% ��������
function [Hmin, HminTime, indHminTime] = calcH(kedTime, ked, sysState, mu)
	sysASize = size(sysState,2);
	lenKedTime = length(kedTime);
	% ��������� � ������ ���
	Hmin = Inf;
	for j = 1 : lenKedTime
		H = 0;
		for k = 1 : sysASize % ���������, ��� ������� � �������� �������� ���������
			% �.�. ������ ������ ��������� ������ = 1, �� ���� ����� ���� ���������� ��������� � ���
			 H = H + mu(k)*( ( sysState(k) - ked(j,k) )^2 );
%			 H = sqrt(H);
		end;
		if Hmin > H
			HminTime = kedTime(j);
			indHminTime = j;
			Hmin = H;
			if lenKedTime == j % HminTime = kedTime(end)
				  break%, exit from function
			end
		end;
	end;
    
% % ��� ��������� �� ���������
% 	Hmin = mu*((sysState - ked').^2)';
%     HminTime = sum((sysState - ked').^2);
%     indHminTime = 1;

%}
