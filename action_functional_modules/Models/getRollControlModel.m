function [sys_cl, sys_eigs] = getRollControlModel(sys)
% x = .001;
% y = .001;
x = .1;
y = .1;
Q = [x 0;
     0 y];
R = 1; % 50
[K, ~, sys_eigs]= lqr(sys.A,sys.B,Q,R);

Ac = [(sys.A-sys.B*K)];
Bc = [sys.B]; 
%?????? т.е. замкнули систему, а вход у неё остался точно такой же????
% да, такой же. Это вход управления, а не возмущения. Референсный сигнал r
% входит с таким же "В": dot(x) = A*x + B*r - B*K*x = (A-B*K)*x + B*r
% В данной системе вход и есть само возмущение. 
Cc = [sys.C];
Dc = [sys.D];

sys_cl = ss(Ac,Bc,Cc,Dc);
sys_cl.StateName = {'omega'; 'domega'};
sys_cl.InputName = {'waveSlope'};
sys_cl.OutputName = {'omega'; 'domega'};

% sys_eigs = eig(sys_cl.A);

% clear x y Q R K Ac Bc Cc Dc

%%
%{
T=0:0.01:10;
U=0.2*ones(size(T));
[Y,T,X]=lsim(sys_cl,U,T);

figure
plot(T,Y)
%}
end