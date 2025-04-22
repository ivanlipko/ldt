function sys = getRollModel(courseAngle)
%% ss модель Бортовой качки
% courseAngle - угол в градусах

ship
shipRoll

% courseAngle = 0; % курсовой угол в градусах
% p.161 Vagushenko
HI_teta = HI_tetaT * sind(courseAngle);
% 1 - (2 * pi * T * 3.62)/(0.715 * 1.62)
%(5.43) Vagushenko
b_teta = omega_teta^2 * HI_teta ;

% модель очень сильно зависит от длины волны, т.к. HI_tetaT = func(waveLength)
% omega'' + 2*omega_teta*Mu_teta * omega' + omega_teta*omega_teta * omega = b_teta * waveSlope

A = [0  1;
     -omega_teta*omega_teta    -2*omega_teta*Mu_teta]; 
B = [0; b_teta];
C = [1 0;
     0 1];
D = [0;0];

sys = ss(A,B,C,D);
sys.StateName = {'omega'; 'domega'};
sys.InputName = {'waveSlope'};
sys.OutputName = {'omega'; 'domega'};

% clear A B C D G

%%
% time = 0: 0.01: 30;
% figure
% step(sys, time)
% [y_impulse, t_impulse] = impulse(sys, time);
% impulse(sys, time);

% p = eig(A)
end