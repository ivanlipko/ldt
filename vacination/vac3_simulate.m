%% Симуляция системы, чтобы посомтреть как она себя ведёт, 
% что представляет из себя

clear
clc
% close all



%% Добавляю пути к скриптам
[path_main, ~, ~] = fileparts(mfilename('fullpath'));
cd(path_main)
addpath(genpath( './../../action_functional_modules' ))

% if 1 
%     return 
% end

%% Модель SERIRS 
[rightSERIRS, getSERIRSinfo, rightSERIRSLinear, rightSERIRSLinearDiseaseFree] = getVacModel3();

Asize = 4;

% Parameters from examples 2.5 and 2.6
alpha = 1/10;
gamma = 1/7;
delta = 1/14;
n = 100;
sigma = 1/7;
omega = 1/90;

beta = 0.4;  % Ro = 2.053
% beta = 0.2;  % Ro = 1.027
% beta = 0.195;  % Ro = 0.5133

[p, eps, Ro, is] = getSERIRSinfo(alpha, beta, gamma, delta, sigma, omega, n)

rng(100)
t = 0:0.01:360*2;
w = zeros(size(t,2), Asize);
w = randn(size(t,2), Asize);
% w = 0.05 * w;

% s e i r
S0 = 90;  % 20, 30, 40, 50, 60, .. 90
% x0 = [S0, (n-S0)/2, (n-S0)/2, 0]
x0 = [S0, (n-S0), 0, 0]

if sum(x0) ~= n
    error("Initial condition irrelevent. Must be sum(x0(1:3)) == n")
end

param.alpha = alpha;
param.beta = beta;
param.sigma = sigma;
param.gamma = gamma;
param.delta = delta;
param.omega = omega;
param.n = n;

out = rkSolve(@rkStepSERIRS, rightSERIRS, Asize, t, x0, param, w);
fprintf('Integration finish\n')

%
figure(1), clf, grid on, hold on, legend('Location', 'south')
subplot(221), grid on, hold on, title('SERIRS Model') %, legend('Location', 'best')
xlabel('s'), ylabel('e')
plot(out(1,:),out(2,:), 'DisplayName', 'Path')
plot(out(1,1),out(2,1), 'o')


subplot(222), grid on, hold on%, legend('Location', 'best')
xlabel('s'), ylabel('i')
plot(out(1,:),out(3,:), 'DisplayName', 'Path')
plot(out(1,1),out(3,1), 'o')

subplot(223), grid on, hold on%, legend('Location', 'best')
xlabel('e'), ylabel('i')
plot(out(2,:),out(3,:), 'DisplayName', 'Path')
plot(out(2,1),out(3,1), 'o')

subplot(224), grid on, hold on%, legend('Location', 'best')
xlabel('i'), ylabel('r')
plot(out(3,:),out(4,:), 'DisplayName', 'Path')
plot(out(3,1),out(4,1), 'o')


figure(2), clf
subplot(211), grid on, hold on, legend('Location', 'best'), title('SERIRS Model')
xlabel('Time, days'), ylabel('')
plot(t,out(1,:), 'DisplayName', 's')
plot(t,out(2,:), 'DisplayName', 'e')

subplot(212), grid on, hold on, legend('Location', 'best')
xlabel('Time, days'), ylabel('')
plot(t,out(3,:), 'DisplayName', 'i')


%%
% 
% x state vector = (s,e,i)

param.alpha = alpha;
param.beta = beta;
param.sigma = sigma;
param.gamma = gamma;
param.delta = delta;
param.omega = omega;
param.n = n;
param.phi = 'none';
param.psi = eps;
param.ro = Ro;

out_linear = rkSolve(@rkStepSVERIRS, rightSERIRSLinear, Asize-1, t, x0([1,2,3]), param, w(:,[1,2,3]));
out_linear(4,:) = n - sum(out_linear([1,2,3],:), 1);
% out_linear = out_linear+p; 
% out_linear([1,2,3],:) = out_linear([1,2,3],:) + p([1,2,3]); 
fprintf('Integration finish\n')

figure(3),  clf
subplot(221), grid on, hold on, title('SERIRS Linear Model')  %, legend('Location', 'best')
xlabel('s'), ylabel('e')
plot(out_linear(1,:),out_linear(2,:), 'DisplayName', 'Path')
plot(out_linear(1,1),out_linear(2,1), 'o')

subplot(222), grid on, hold on%, legend('Location', 'best')
xlabel('s'), ylabel('i')
plot(out_linear(1,:),out_linear(3,:), 'DisplayName', 'Path')
plot(out_linear(1,1),out_linear(3,1), 'o')

subplot(223), grid on, hold on%, legend('Location', 'best')
xlabel('e'), ylabel('i')
plot(out_linear(2,:),out_linear(3,:), 'DisplayName', 'Path')
plot(out_linear(2,1),out_linear(3,1), 'o')

subplot(224), grid on, hold on%, legend('Location', 'best')
xlabel('i'), ylabel('r')
plot(out_linear(3,:),out_linear(4,:), 'DisplayName', 'Path')
plot(out_linear(3,1),out_linear(4,1), 'o')

figure(4), clf
subplot(211), grid on, hold on, legend('Location', 'best'), title('SERIRS Linear Model')
xlabel('Time, days'), ylabel('')
plot(t,out_linear(1,:), 'DisplayName', 's') 
plot(t,out_linear(2,:), 'DisplayName', 'e')

subplot(212), grid on, hold on, legend('Location', 'best')
xlabel('Time, days'), ylabel('')
plot(t,out_linear(3,:), 'DisplayName', 'i')



%% SVERIRS model
clc
[~, ~, ~, ~, rightSVERIRS, getSVERIRSinfo, rightSVERIRSLinear] = getVacModel3();
Asize = 5;

% Parameters from examples 3.2
alpha = 1/10;
% beta = 1/5;   % example 3.2
beta = 9/10;   % example 3.3. Endemic equilibirum in p = (s,e,i,r,v) = (21,3,3,66,7)
gamma = 1/7;
delta = 1/14;
n = 100;
sigma = 1/7;
omega = 1/90;
phi = 1/360;
psi = 1/180;
ro = 1/10;

[p, eps, Ro,is] = getSVERIRSinfo(alpha, beta, gamma, delta, sigma, omega, n, psi, phi, ro)

t = 0:0.01:36;
w = zeros(size(t,2), Asize);
w = randn(size(t,2), Asize);
% % w = 0.05 * w;

% s e i r v
S0 = 90;  % 20, 30, 40, 50, 60, .. 90
% x0 = [S0, (n-S0)/2, (n-S0)/2, 0, 0]
% x0 = [S0, (n-S0), 0, 0, 0]
% x0 = [n, 0, 0, 0, 0]
x0 = [21,3,3,66,7];

if sum(x0) ~= n
    error("Initial condition irrelevent. Must be sum(x0(1:3)) == n")
end

param.alpha = alpha;
param.beta = beta;
param.sigma = sigma;
param.gamma = gamma;
param.delta = delta;
param.omega = omega;
param.n = n;
param.phi = phi;
param.psi = psi;
param.ro = ro;

out = rkSolve(@rkStepSVERIRS, rightSVERIRS, Asize, t, x0, param, w);
fprintf('Integration finish\n')

%
figure(5), clf
subplot(221), grid on, hold on, title('SVERIRS Model') %, legend('Location', 'best')
xlabel('s'), ylabel('e')
plot(out(1,:),out(2,:), 'DisplayName', 'Path')
plot(out(1,1),out(2,1), 'o')

subplot(222), grid on, hold on%, legend('Location', 'best')
xlabel('s'), ylabel('i')
plot(out(1,:),out(3,:), 'DisplayName', 'Path')
plot(out(1,1),out(3,1), 'o')

subplot(223), grid on, hold on%, legend('Location', 'best')
xlabel('e'), ylabel('i')
plot(out(2,:),out(3,:), 'DisplayName', 'Path')
plot(out(2,1),out(3,1), 'o')

subplot(224), grid on, hold on%, legend('Location', 'best')
xlabel('i'), ylabel('r')
plot(out(3,:),out(4,:), 'DisplayName', 'Path')
plot(out(3,1),out(4,1), 'o')


figure(6), clf
subplot(211), grid on, hold on, legend('Location', 'best'), title('SVERIRS Model')
xlabel('Time, days'), ylabel('')
plot(t,out(1,:), 'DisplayName', 's')
plot(t,out(2,:), 'DisplayName', 'e')
plot(t,out(3,:), 'DisplayName', 'i')

subplot(212), grid on, hold on, legend('Location', 'best')
xlabel('Time, days'), ylabel('')
plot(t,out(4,:), 'DisplayName', 'r')
plot(t,out(5,:), 'DisplayName', 'v')


%%
out_linear = rkSolve(@rkStepSVERIRS, rightSVERIRSLinear, Asize-1, t, x0([1,2,3,5]), param, w(:,[1,2,3,5]));

figure(7),  clf
subplot(221), grid on, hold on%, legend('Location', 'best')
xlabel('s'), ylabel('e')
plot(out_linear(1,:),out_linear(2,:), 'DisplayName', 'Path')
plot(out_linear(1,1),out_linear(2,1), 'o')

subplot(222), grid on, hold on%, legend('Location', 'best')
xlabel('s'), ylabel('i')
plot(out_linear(1,:),out_linear(3,:), 'DisplayName', 'Path')
plot(out_linear(1,1),out_linear(3,1), 'o')

subplot(223), grid on, hold on%, legend('Location', 'best')
xlabel('e'), ylabel('i')
plot(out_linear(2,:),out_linear(3,:), 'DisplayName', 'Path')
plot(out_linear(2,1),out_linear(3,1), 'o')

subplot(224), grid on, hold on%, legend('Location', 'best')
xlabel('i'), ylabel('r')
plot(out_linear(3,:),out_linear(4,:), 'DisplayName', 'Path')
plot(out_linear(3,1),out_linear(4,1), 'o')


%%

%%
function out = rkSolve(rkStep, rfun_handle, Asize, t, x0, param, w)
xvec = x0' ; % Вектор состояния. Начальные условия
out = zeros(Asize,length(t));  % История вектора состояния
out(:,1) = xvec;
dh = t(2)-t(1);

% переменные симулирования
% pbCount = length(t) / 10; % прогрессБар - всего шагов
% pbStep = 0;                     % прогрессБар первый шаг
% pb = 10;                        % прогрессБар - текущий шаг
k = 2;
% tic
for tk = t(2:end)
%     pbStep = pbStep + 1;
%     if pbStep > pbCount
%         fprintf('%1.0f..', pb);
%         pbStep = 0;
%         pb = pb + 10;
%     end
%     
    xvec = rkStep(tk, dh, xvec, rfun_handle, param, w(k,:));
    out(:,k) = xvec; % when vector subscript to (:,)

    k = k + 1;
end
% fprintf('%1.0f.. The end. ', pb);
% toc
end

%%
function X = rkStepSERIRS(t,h,X, right, param, wk)
    alpha = param.alpha;
    beta = param.beta;
    gamma = param.gamma;
    omega = param.omega;
    sigma = param.sigma;
    delta = param.delta; 
    n = param.n;
    
	h2=0.5*h;
	h6=0.166666666*h;

    F=right(t,X, wk, alpha, beta, gamma, delta, omega, sigma, n);
    Xr=X+h2*F;
    s=F;    
    F=right(t+h2,Xr, wk, alpha, beta, gamma, delta, omega, sigma, n);
    
	Xr=X+h2*F;
    s=s+2*F;	
    F=right(t+h2,Xr, wk, alpha, beta, gamma, delta, omega, sigma, n);

    Xr=X+h*F;
	s=s+2*F;

    F=right(t+h2,Xr, wk, alpha, beta, gamma, delta, omega, sigma, n);
    
	X=X+h6*(s+F);
end


%%
function X = rkStepSVERIRS(t,h,X, right, param, wk)
    alpha = param.alpha;
    beta = param.beta;
    gamma = param.gamma;
    omega = param.omega;
    sigma = param.sigma;
    delta = param.delta; 
    n = param.n;
    phi = param.phi;
    psi = param.psi;
    ro = param.ro;
    
	h2=0.5*h;
	h6=0.166666666*h;

    F=right(t,X, wk, alpha, beta, gamma, delta, omega, sigma, n, phi, psi, ro);
    Xr=X+h2*F;
    s=F;    
    F=right(t+h2,Xr, wk, alpha, beta, gamma, delta, omega, sigma, n, phi, psi, ro);
    
	Xr=X+h2*F;
    s=s+2*F;	
    F=right(t+h2,Xr, wk, alpha, beta, gamma, delta, omega, sigma, n, phi, psi, ro);

    Xr=X+h*F;
	s=s+2*F;

    F=right(t+h2,Xr, wk, alpha, beta, gamma, delta, omega, sigma, n, phi, psi, ro);
    
	X=X+h6*(s+F);
end
