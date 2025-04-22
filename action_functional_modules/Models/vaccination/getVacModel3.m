%% Модель Вакцинации из статьи
% TWO NEW COMPARTMENTAL EPIDEMIOLOGICAL MODELS ANDTHEIR EQUILIBRIA by 
% JONAS BALISACAN, MONIQUE CHYBA, AND COREY SHANBROM
% 
% https://arxiv.org/pdf/2109.01738.pdf
% 
% Time is day
% 
% Parameters from example 2.5
% alpha = 1/10;
% gamma = 1/7;
% delta = 1/14;
% n = 100;
% sigma = 1/7;
% omega = 1/90;
% 
% beta = 0.4;  % Ro = 2.053
% 
%  Endemic eqilibrium at point
% p = (s,e,i,r) = n/Ro (1, omega * eps, sigma*omega * eps / gamma, (sigma+delta) * eps )
% where eps = 1/(sigma+delta) * (beta * (alpha*gamma+sigma) - gamma * (delta+sigma)) / (sigma*omega+gamma*(delta+sigma+omega))
% Endemic equilibrium is realistic if condition (9) solved:
%  beta * (alpha*gamma+sigma) > gamma*(sigma+delta)
% which is eqivalent to Ro > 1
% 
% 
% ЭНДЕМИЧНЫЙ. (endemic) - часто возникающий в каком-либо определенном регионе 
% или у какой-либо определенной части населения: данный термин применяется 
% по отношению к заболеваниям, которые обычно или постоянно присутствуют в 
% данном сообществе людей.
% 

function [rightSERIRS_, getSERIRSinfo_, rightSERIRSLinearEndemic_, rightSERIRSLinearDiseaseFree_,...
    rightSVERIRS_, getSVERIRSinfo_, rightSVERIRSLinear_, ...
    getSysSERIRSLinearEndemic_, getSysSERIRSLinearDiseaseFree_, getSysSVERIRSLinDiseaseFree_, ...
    relocateColumnsinSys_, ...
    rkSolve_, rkStepSERIRS_, rkStepSVERIRS_] = getVacModel3()
%%
rkSolve_ = @rkSolve;
rkStepSERIRS_ = @rkStepSERIRS;
rkStepSVERIRS_ = @rkStepSVERIRS;

rightSERIRS_ = @rightSERIRS;
getSERIRSinfo_ = @getSERIRSinfo;
rightSERIRSLinearEndemic_ = @rightSERIRSLinearEndemic;
rightSERIRSLinearDiseaseFree_ = @rightSERIRSLinearDiseaseFree;

getSysSERIRSLinearEndemic_ = @getSysSERIRSLinearEndemic;
getSysSERIRSLinearDiseaseFree_ = @getSysSERIRSLinearDiseaseFree;

rightSVERIRS_ = @rightSVERIRS;
getSVERIRSinfo_ = @getSVERIRSinfo;
rightSVERIRSLinear_ = @rightSVERIRSLinear;

getSysSVERIRSLinDiseaseFree_ = @getSysSVERIRSLinDiseaseFree;

relocateColumnsinSys_ = @relocateColumnsinSys;
end


%% SE(R)IRS model
% Formulas (1)-(4)
function dx = rightSERIRS(t,x, w, alpha, beta, gamma, delta, omega, sigma, n)
    s = x(1);
    e = x(2);
    i = x(3);
    r = x(4);
    
    ds = -beta * s * (i+alpha*e)/n + omega * r + w(1);
    de =  beta * s * (i+alpha*e)/n - (sigma + delta) * e + w(2);
    di = sigma*e - gamma*i + w(3);
    dr = delta*e + gamma*i - omega*r + w(4);
    
    dx = [ds; de; di; dr];
end


% Formula for M из раздела 2.1
% Lineariation at point p
function dx = rightSERIRSLinearEndemic(t,x, w, alpha, beta, gamma, delta, omega, sigma, n, phi, psi, ro)
% Здесь ro это Ro
eps = psi;
    M = [-eps*omega*(delta+sigma)-omega,   -alpha*beta/ro-omega,    -beta/ro-omega;
              eps*omega*(delta+sigma),     -beta*sigma/(gamma*ro),   beta/ro;
              0,                                  sigma,            -gamma];

    dx = M*x + w';
end

% 
% (s,e,i,r)
function sys = getSysSERIRSLinearEndemic(alpha, beta, gamma, delta, omega, sigma, eps, Ro, B,C,D)
if nargin <= 8
    B = [-1; 0;1]; 
    C = eye(3);
    D = [];
end
    M = [-eps*omega*(delta+sigma)-omega,   -alpha*beta/Ro-omega,    -beta/Ro-omega;
              eps*omega*(delta+sigma),     -beta*sigma/(gamma*Ro),   beta/Ro;
              0,                                  sigma,            -gamma];
    sys = ss(M,B,C,D);
    sys.StateName = {'S'; 'E'; 'I'};
    sys.StateUnit = {'Восприимчивые, чел'; 'Инфецированные без симптомов, чел'; 'Инфецированные с симптомами, чел'};
%     syslin.StateName = {'Восприимчивые'; 'Инфецированные без симптомов'; 'Инфецированные с симптомами'; 'Выздоровевшие'};
end


% Formula for N из раздела 2.2
% Disease-free equilibrium at (s,e,i,r) = (n,0,0,0)
function dx = rightSERIRSLinearDiseaseFree(t,x, w, alpha, beta, gamma, delta, omega, sigma, n, phi, psi, ro)
    A = [-omega,   -alpha*beta-omega,   -beta-omega;
              0,     alpha*beta-delta-sigma,  beta;
              0,  sigma,  -gamma;];

    dx = A*x + w';
end

% Analysis of disease-free equilibria.
% disease-free equilibrium at (S, E, I, R) = (n, 0, 0, 0)
function sys = getSysSERIRSLinearDiseaseFree(alpha, beta, gamma, delta, omega, sigma, B,C,D)
if nargin < 7
    B = [-1; 0;1]; 
    C = eye(3);
    D = [];
end
    N = [-omega,   -alpha*beta-omega,   -beta-omega;
              0,     alpha*beta-delta-sigma,  beta;
              0,  sigma,  -gamma;];
    sys = ss(N,B,C,D);
    sys.StateName = {'S'; 'E'; 'I'};
    sys.StateUnit = {'Восприимчивые, чел'; 'Инфецированные без симптомов, чел'; 'Инфецированные с симптомами, чел'};
%     syslin.StateName = {'Восприимчивые'; 'Инфецированные без симптомов'; 'Инфецированные с симптомами'; 'Выздоровевшие'};
end


function [p, eps, Ro,is] = getSERIRSinfo(alpha, beta, gamma, delta, sigma, omega, n)
% p = (s,e,i,r)  -- eqilibrium
    eps = 1/(sigma+delta) * (beta * (alpha*gamma+sigma) - gamma * (delta+sigma)) / (sigma*omega+gamma*(delta+sigma+omega));
    Ro = (alpha*gamma+sigma)/(delta+sigma) * beta/gamma;
    p = n/Ro * [1; omega * eps; sigma*omega * eps / gamma; (sigma+delta) * eps];
    is = beta * (alpha*gamma+sigma) > gamma*(sigma+delta);
end


%%
function sys = relocateColumnsinSys(sys, cols)
    sys.A = sys.A(cols,cols);
%     sys.B = sys.B(cols,cols);
    sys.B = sys.B(cols);
    sys.C = sys.C(cols,cols);
%     sys.D = sys.D(cols,cols);
    sys.D = sys.D(cols);
    sys.StateName = sys.StateName(cols);
    sys.StateUnit = sys.StateUnit(cols);
end


%% SVE(R)IRS model
% Formulas (14)-(18)
function dx = rightSVERIRS(t,x, w, alpha, beta, gamma, delta, omega, sigma, n, phi, psi, ro)
    s = x(1);
    e = x(2);
    i = x(3);
    r = x(4);
    v = x(5);
    
    ds = -beta * s * (i+alpha*e)/n + omega * r - phi * s + psi * v + w(1);
    de =  beta * s * (i+alpha*e)/n - (sigma + delta) * e + ro*beta*v*(i+alpha*e)/n + w(2);
    di = sigma*e - gamma*i + w(3);
    dr = delta*e + gamma*i - omega*r + w(4);
    dv = -ro*beta*v*(i+alpha*e)/n + phi*s - psi*v + w(5);
    
    dx = [ds; de; di; dr; dv];
end

% Formula for N в разделе 3.2. Как-будто что-то не то с моделью
function dx = rightSVERIRSLinear(t,x, w, alpha, beta, gamma, delta, omega, sigma, n, phi, psi, ro)
    A = [-phi-omega,   alpha*beta*psi/(phi+psi)-omega,   -beta*psi/(phi+psi)-omega,   psi-phi;
              0,     -delta-sigma-alpha*beta*(psi+ro*psi)/(phi+psi),  beta*(psi+ro*phi)/(phi+psi),  0;
              0,  sigma,  -gamma,    0;
              phi,  -alpha*beta*ro*phi/(phi+psi),   -beta*ro*phi/(phi+psi),   -psi];

          % пробовал искать ошибку, но тут её нет:
          % 0,     -delta-sigma+alpha*beta*(psi+ro*psi)/(phi+psi),  beta*(psi+ro*phi)/(phi+psi),  0;
    dx = A*x + w';
end

function [p1, Ro, is, eps] = getSVERIRSinfo(alpha, beta, gamma, delta, sigma, omega, n, psi, phi, ro)
% p = (s,e,i,r,v)  -- disease-free eqilibrium
eps='none';
p1 = [ psi*n/(phi+psi), 0, 0, 0, phi*n/(phi+psi)];
is = 'none';
    Ro = (alpha*gamma+sigma)/(delta+sigma) * (psi + ro*phi)/(psi+phi) * beta/gamma;
    is = beta * (alpha*gamma+sigma) > gamma*(sigma+delta);
end

% Disease-free linear system . Chapter 3.2
% disease-free equilibrium at (S, E, I, R, V)
% State (S, E, I, V)
function sys = getSysSVERIRSLinDiseaseFree(alpha, beta, gamma, delta, omega, sigma, psi, phi, ro, B,C,D)
if nargin < 10
    B = [-1; 0; 0; 1];   % вакцинируем
    C = eye(4);
    D = [];
end
    N = [-phi-omega,   alpha*beta*psi/(phi+psi)-omega,   -beta/(phi+psi)-omega,  psi - omega;
              0,     -delta - sigma - alpha*beta*(psi+ro*psi)/(phi+psi),  beta*(psi+ro*psi)/(phi+psi), 0;
              0,  sigma,  -gamma  0;
              phi, -alpha*beta*ro*phi/(phi+psi),  -beta*ro*phi/(phi+psi),   -psi];
    sys = ss(N,B,C,D);
    sys.StateName = {'S'; 'E'; 'I'; 'V'};
    sys.StateUnit = {'Восприимчивые, чел'; ...
        'Инфецированные без симптомов, чел'; ...
        'Инфецированные с симптомами, чел'; ...
        'Вакцинировано, чел'};
end





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


