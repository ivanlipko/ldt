%% Модель вакцинации
% из статьи Dani Suandi, Optimal Control Problem of Vacination for The
% Spread of Measles Diseases Model. 2018
% 
% Measles - Корь, краснуха 
% 
% s - susceptible group (восприимчивая группа)
% i - infected group (инфицированные)
% r - recovery group (выздоравлившие, повторно заболеть не могут)
% 
% According to Arcticle Ro < 1
% mu = 0.01;  % natural mortality, time^-1
% beta = 0.2;  % infection rate, time^-1
% gamma = 0.03;  % natural recovery rate, time^-1
% u = 0.06;  % vaccination rate

% According to Arcticle Ro > 1
% mu = 0.01;  % natural mortality, time^-1
% beta = 0.3;  % infection rate, time^-1
% gamma = 0.01;  % natural recovery rate, time^-1
% u = 0.05;  % vaccination rate


function [fright, fright_WithoutN] = getVacModel1()
% особые точки
% ps = [gm/dl al/bt];

fright = @right;
       
fright_WithoutN = @right_WithoutN;

end
%%

% In this system u \in [0;1) is vaccination control rate
function dx = right(t,x, u, beta, mu, gamma)
    s = x(1);
    i = x(2);
    r = x(3);
    N = x(4);
    lambda = mu * N;
    ds = lambda - beta*(1-u) * s*i/N - (u+mu) * s;
    di = beta*(1-u) * s*i/N - (gamma + mu) * i;
    dr = u * s + gamma * i - mu * r;
    N = s + i + r;
    
    dx = [ds; di; dr; N];
end

%%
%  formula (2) from the article
function dx = right_WithoutN(t,x, u, beta, mu, gamma)
    s = x(1);
    i = x(2);

    ds = mu - beta*(1-u) * s*i - (u+mu) * s;
    di = beta*(1-u) * s*i - (gamma + mu) * i;
    
    dx = [ds; di];
end

