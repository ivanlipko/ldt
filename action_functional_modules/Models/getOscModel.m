%% Модель осцилятора второго порядка.
%  
function sys = getOscModel()

% A  = [0 1;
%      -1 0];
% G = [ 0;  1];


a1 = .5;
a2 = .05;
A  = [0         1         
     -a1      -a2];

g1 = .01;
g2 = 1;
G = [ g1  0     
      0   g2];

% модифицировано: упрощено, один вход
% g1 = .00001;
% G = [ g1;  
%     g1];
  
C = eye(size(A));
D = zeros(size(G));


% D = [ 0;  0];

sys = ss(A,G,C,D);

sys.Name = 'Linear Oscilation Model (from Dubovik example)';
sys.StateName = {'One', 'Two' };
% sys.InputName = {'In1', 'In2'};
sys.InputName = {'In1'};
sys.OutputName = sys.StateName;
