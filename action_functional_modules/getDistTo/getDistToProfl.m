% Оценка близоститекущего состояния системы к Профилю уклонения
%
%  curState - вектор состояния объекта в текущий момент времени (!!!)
%  mu - вектор векторов весов
%     size(mu) == size(sysState)
%
%  Обновления:
%   2022.05.31
%       Переименовал indMinDist в minDistInd
%   2022.04.28
%       Проверка размерности времени и траектории path
%   2020.08.28
%       ДОбавил в Выходной параметр структуру, содержащую теже параметры.
%       Рефакторинг
%   2018.10.27
%       Рефакторинг. Переименование модуля из calcH в getDistToProfl
%   2018.04.18
%       Изменение. Учитывается размер системы и вектора mu: он изменяется в
%       зависимости от их размеров. ones(size(kedTime,2), size(mu,2))
%   2018.03.20
%       Ускорение работы скрипта в 31 раз за счёт использования векторных
%       операций. Формула таже.

function [minDist, minDistTime, minDistInd, v] = getDistToProfl(time, path, cstate, mu)
assert( sum(size(mu) == size(path(1,:)) ) == 2, 'sum(size(mu) == size(prfl(1,:))) == 2 ' )
assert( size(cstate,2) == size(path,2) , 'must be TRUE size(cstate,2) == size(path,2)' )
assert( size(time,2) == size(path,1) , 'must be TRUE size(time,2) == size(path,1) ' )

dist = mu .* ones(size(time,2), size(mu,2)) .* ( (cstate - path).^2 );   % euclid
% dist = mu .* ones(size(time,2), size(mu,2)) .* sum( abs(cstate - path),2 );     % manhatten
% sum(abs(x-y))

% sqrt??? тогда будет как норма Эвклида
[minDist, minDistInd ] = min(sum(dist, 2));
minDistTime = time(minDistInd);

v.minDist = minDist;
v.minDistTime = minDistTime;
v.indMinDist = minDistInd;
end



%{
% оригинал
function [Hmin, HminTime, indHminTime] = calcH(kedTime, ked, sysState, mu)
	sysASize = size(sysState,2);
	lenKedTime = length(kedTime);
	% сравнение с точкой КЭД
	Hmin = Inf;
	for j = 1 : lenKedTime
		H = 0;
		for k = 1 : sysASize % учитывать, что работаю с вектором векторов состояния
			% т.к. сейчас вектор состояния сстемы = 1, то ищем среди всех траекторий ближайщую к КЭД
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
    
% % для сравнения по тректории
% 	Hmin = mu*((sysState - ked').^2)';
%     HminTime = sum((sysState - ked').^2);
%     indHminTime = 1;

%}
