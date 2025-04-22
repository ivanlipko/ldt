% Модель катамарана 
% length 90m, beam 25.96m, wareline Length 85.655m
% dispalcement 740t
% draft 2.6m
%
% Пример использования:
%   angleDegPeterb = 90;  % лагом к волне
%   %angleDegPeterb = 0;  % носом к волне
%   cols = [6, 5, 3, 1, 4, 2];  % roll first 
%   [sys, sys_feedbacked] = getKatamaran(cols, angleDegPeterb);
% 
% В статье условия симулирования:
% catamaran sail on 40 knots
% level sea condition 5
%  wave height 3.25
% wave approaching angle 180 degrees
% 
% dot(x) = A*x+ B*u + G*w
% 'u' is control input
% 'w' is environment disturbance
% 
% Сделать:
%   перестановка входов в зависимости от cols
% Обновления
%   27.10.2018 Вложил внутрь LQ-матрицы. Матрица внешних возмущений G 
%       зависит от курсовго угла на генерального направление волны. Рефаторинг
%   19.10.2018 Изменения динамики; рефакторинг
%   02.10.2018 Изменения кода перестановки названий состояний
%   23.04.2018 Перестановка столбцов матрицы Q для LQ-регулятора
%   19.04.2018 Перестановка столбцов матриц управления и объекта


function [katamaran, katamaranFb] = getKatamaran(cols, angleDegPeterb)
% из статьи H? Output Feedback Control Method for the 
% T-Foil of the Wave Piercing Catamaran by Songtao Zhang.   (07053443.pdf)
% В статье при моделировании указано, что волна идёт под углом 180 градусов (в корму)
% A = [-0.9282 -25.5640 -11.1113 -23.4295;
%     0.0514 -0.5231 0.2503 -12.1951;
%     1 0 0 0;
%     0 1 0 0];
% B = [0.6670; 
%     0.0411;
%     0;
%     0];

% Из статьи Application of Model Predictive Control Technique for Wave 
% Piercing Catamarans Ride Control System by Lihua Liang, Jia Yuan, and 
% Songtao Zhang (07558652.pdf). Дополнительно можно посмотреть:
% Ride Control Method of Wave-Piercing Catamaran with T-foil and Flaps by
% Songtao Zhang, Shuanglin Li, Lihua Liang, Mingxiao Sun (06885750.pdf)

% State Names are Heave Vel; Pitch Vel; Heave; Pitch;
% stNames = {'Heave Vel'; 'Pitch Vel'; 'Heave'; 'Pitch'};
% stUnits = {'m/s'; 'deg/s'; 'm'; 'deg'};
% A = [-0.9073 -25.1097 -14.1503 -17.4945;
%     0.0514  -0.503  0.2442  -12.4;
%     1   0   0   0;
%     0   1   0   0    ];
% B = [0.0082 0.0000083;
%     -0.00016    0.000017;
%     0   0;
%     0   0    ];
% C = [1 0 0 0;
%     0 1 0 0;
%     0 0 1 0;
%     0 0 0 1]; 
% D = zeros(4,4);
% 
% cols = [6, 5, 3, 1, 4, 2]; % roll first

% stNames = {'Heave_Vel'; 'Pitch_Vel'; 'Heave'; 'Pitch'; 'Roll_Vel'; 'Roll'};
% stUnits = {'m/s'; 'deg/s'; 'm'; 'deg'; 'deg/s'; 'deg'};

stNames = {'$$\dot\eta$$'; '$$\dot\zeta$$'; '$$\eta$$'; '$$\zeta$$'; '$$\dot\theta$$'; '$$\theta$$'};
% stNames = {'d\eta / dt'; 'd\zeta / dt'; '\eta'; '\zeta'; 'd\theta / dt'; '\theta'};
% stNames = {'Heave_Vel'; 'Pitch_Vel'; 'Heave'; 'Pitch'; 'd\theta / dt'; '\theta'};
stUnits = {'Heave Velocity, m/s'; 'Pitch Velocity, deg/s'; 'Heave, m'; 'Pitch, deg'; 'Roll Velocity, deg/s'; 'Roll, deg'};

% stNames = {'Скорость вертикальной качки'; 'Скорость дифферента'; 'Вертикальная качка'; 'Дифферент'; 'Скорость бортовой качки'; 'Бортовая качка'};
% stUnits = {'м/с'; '\circ/с'; 'м'; '\circ'; '\circ/с'; '\circ'};

A = [-0.9073 -25.1097 -14.1503 -17.4945   10^-3 0;
    0.0514   -0.503    0.2442  -12.4      10^-3 0;
    1        0         0        0         0    0 ;
    0        1         0        0         0    0 ;
    10^-3    10^-3     -10^-2    -10^-2    -5 -15;
    0        0         0        0         1    0];
% A = [-0.9073 -25.1097 -14.1503 -17.4945   10^-3 0;
%     0.0514   -0.503    0.2442  -12.4      10^-3 0;
%     1        0         0        0         0    0 ;
%     0        1         0        0         0    0 ;
%    0 0 0 0     -.2 -.5;
%     0        0         0        0         1    0];
B = [0.0082    0.0000083;
    -0.00016    0.000017;
    0   0;
    0   0;
    0.002 0.0000083;
    0   0];

% G = B; % Потому что урпавление и возмущение поступают в качестве сил - 
% одна и та же размерность. Но возмущения зависят от курсового угла на
% генеральный фронт волны
G = [0.0082*cosd(angleDegPeterb)     0.0000083*cosd(angleDegPeterb);
    -0.00016*cosd(angleDegPeterb)    0.000017*cosd(angleDegPeterb);
    0   0;
    0   0;
    0.0082*sind(angleDegPeterb) 0.0000083*sind(angleDegPeterb);
    0   0];

BG = [B, G];

%% Перестановки
stNames = stNames(cols);
stUnits = stUnits(cols);
A = A(cols, cols);
BG = BG(cols, :);
B = B(cols,:);
G = G(cols,:);
% перестанвока входов.

%%
% katamaran = ss(A, BG, eye(6,6), zeros(6,4));   % почему BG??
katamaran = ss(A, B, eye(6,6), zeros(6,2));
katamaran.Name = 'Catamaran';
katamaran.StateName = stNames;
katamaran.StateUnit = stUnits;
% katamaran.InputName = {'F_control', 'M_control', 'F_disturbance', 'M_disturbance'};
katamaran.InputName = {'Force', 'Moment'};
katamaran.OutputName = katamaran.StateName;

%% LQ-регулятор
% Строим управление для судна (обычно расположенного лагом к волне).
% При повороте судна обычно не меняют его коэффициенты в ОС, поэтому
% логично не менять всю модель катамарана, а менять коэффициенты для
% действующих возмущений.
Q = [10^1 0 0 0   0 0;
    0 10^1 0 0   0 0 ;
    0 0 4*10^6 0 0 0 ;
    0 0 0 10^1   0 0;
    0 0 0 0      10^4 0;
    0 0 0 0      0  10^4];
Q = 10 * Q;  %  чделал слабое управление, плтомк что слишком силтно  давит качку.  2022.05.20
Q = Q(cols, cols);
R = [0.01 0 ;
    0 0.01];
N = zeros(6,2);

K = lqr(A,B, Q, R, N);
Afb = A-B*K;

katamaranFb = ss(Afb, G, eye(6,6), zeros(6,2));
katamaranFb.Name = 'Catamaran with LQ feedback';
katamaranFb.StateName = stNames;
katamaranFb.StateUnit = stUnits;
katamaranFb.InputName = {'Force', 'Moment'} ;  % {'F disturbance', 'M disturbance'};
katamaranFb.OutputName = katamaranFb.StateName;

end
