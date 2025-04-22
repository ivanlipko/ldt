%% Модель системы первого порядка
%  
function sys = get1Dsys()
a = 2;
b = .0001;

C = 1;
D = 0;

sys = ss(a,b,C,D);

% sys.Name = 'Linear Oscilation Model (from Dubovik example)';
% sys.StateName = {'One', 'Two' };
% sys.InputName = {'In1', 'In2'};
% sys.InputName = {'In1'};
% sys.OutputName = sys.StateName;
