%% Модуль получения профиля и соответсвующих значений функционала действия
%
% Изменения:
% 19.08.2020
%     Переход на другую версию расчёта
% 05.09.2019
%     Рефакторинг pfl.af теперь pfl.aaf . Потому что это А-профиль
% 15.06.2019
%     Изменил выход функции.
%     Результат работы фукнции это структура, которая содержит время,
%     таректории и проч.
% 25.06.2019
%     Рефакторинг
% 26.10.2018
%     Убрал не используемые параметры на выходе, лишние строки комментариев.
% 19.10.2018
%     Рефакторинг
% 02.10.2018
%     Изменил кодировку документа на UTF-8
%     Попытка построить профиль из текущего положения
% 17.04.2018
%     Добавил перестановку столбцов матрицы системы А и матрицы/вектора В
%        управления, если это необходимо.
%     Рефакторинг
%
% % для зеркальных траекторий - лучше если КЭД небудет перескаться

%% TODO
% - Может быть лучше, если на вход будет SS-система sys, время и т.д.,
% тогдаменьше параметров будет передаваться, время будет сразу
% регулироваться снаружи, а внутри использоваться

%%
% function [ssG_out, time, prfle, psiTf, ...
%     prfPeterb, covMatrix, sysA] = ...
%     getProfile(sys, sysG, xTf, kedTf, d_time, cols, ssPeterb)
% function [time, path, peterb, af] = ...
%     getProfile(sys, sysG, xTf, Tf, dT, cols, ssPeterb)
function [pfl] = ...
    getProfile(sys, sysG, xTf, Tf, dT, cols, ssPeterb)

ASize = size(sys.A,1);

%% Определяем ковариационную матрицу
% Для этого строим модель ssModel с учётом или без учёта внешних возмущений
% (по диагонали стоят дисперсии)
% Это нужно для построения такой траектории, которая наименее затратно для
% возмущения выведит объект управления к точке xTf

if nargin == 7
    % Модель с учётом модели возмущений
    ssPeterbA = ssPeterb.A;
    ssPeterbC = ssPeterb.C;
    ssPeterbG = ssPeterb.G;
    sysAPeterbSize = size(ssPeterbA,1);
    sysA = [
        sys.A    sysG*ssPeterbC(1,:)
        zeros(sysAPeterbSize, ASize)   ssPeterbA  ];
    % ssG = [
    %     zeros(sizeofAssFeedbacked(1), 1)
    %     ssPeterbG];
    sysG = [ % ssG_out = ?
        sysG
        ssPeterbG];
    
    % Перестановка строк матрицы управления
    %  надо делать вручную. Пока не автоматизировал
    % Переменная для xTf должна быть в первой строке. Иначе надо будет делать
    %     перестановку строк и столбцов
    %         sysA = sysA([4,1,2,3,5,6,7,8,9,10,11],[4,1,2,3,5,6,7,8,9,10,11]);
    %         ssG_out = sysG([4,1,2,3,5,6,7,8,9,10,11]);
    ASize = ASize + sysAPeterbSize;  % уточняем размер матрицы А
    
    % для проверки размеров матриц
    %  ------------------------------------------------------------------------
    % [ sys.A    sysG*ssPeterbC(1,:) ];
    % [ zeros(sysAPeterbSize, sysASize-1)   ssPeterbA  ]    ;
    % [
    %     sys.A;%(1:4,1:4) ;
    %     zeros(sysAPeterbSize, sysASize )  ];%-1)  ];
    % [
    %     sysG_out*ssPeterbC(1,:)
    %     ssPeterbA  ];
    %  ------------------------------------------------------------------------
else
    % Модель без учёта модели возмущений
    
    if exist('cols', 'var')
        % Перестановка строк матрицы управления
        % Переменная для xTf должна быть в первой строке.
        sysA = sys.A(cols,cols);
        sysG = sysG(cols, :);
    else
        %  без перестановки, 17.04.2018
        sysA = sys.A;
    end
    % State of fullSystemA is
    % [
    %   ...
    % ]
end

covMat = lyap(sysA, sysG*sysG');

%% Вычисление профиля при псиf
Dz = covMat(1,2:ASize);
D1 = covMat(1,1);
xTf = [xTf;  Dz'*xTf/D1];

% psiTf = inv(covMatrix) * xTf;
psiTf = covMat\ xTf;
% SF = psiTf' * covMatrix * psiTf; - эт может быть важным

%% Симулирование. Построение профиля
% симулирование можно вести и за пределами прогнозируемого диапазона - это даст пищу для размышлений
% kedTime = kedTf - TT : 0.5*keddT : kedTf; % аналог -оо..kedTf % оригинальный подход. Ни разу его не исопльзовал. Наврено удалить надо.
time = 0 : dT : Tf;
apath = zeros(ASize, length(time)); % SizeofFullSystemA-1
aaf = zeros( 1, length(time));
for i=1:length(time)
%     pfl(i,:) = expm( sysA * (time(i)-Tf) ) * (xTf - covMat*psiTf) + ...
%         + covMat * expm( sysA' * (Tf-time(i)) ) * psiTf;
    apath(:,i) = covMat * expm( sysA' * (Tf-time(i)) ) * psiTf;
    
    J0 =  covMat*expm((Tf-time(i))*sysA') - expm((time(i)-Tf)*sysA)*covMat; %(45)
    aaf(i) = -.5* psiTf'*expm((Tf-time(i))*sysA)*J0*psiTf; %(57);
end

ptb = zeros(size(sysG,2), length(time)); % SizeofFullSystemA-1
for i=1:length(time)
    ptb(:,i) =  sysG' * expm( sysA' * (Tf-time(i)) ) * psiTf;
end

%%
pfl.time = time;
pfl.apath = apath;
pfl.ptb = ptb;
pfl.aaf = aaf;

