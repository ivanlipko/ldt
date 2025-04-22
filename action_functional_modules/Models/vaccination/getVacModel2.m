%% Модель вакцинации
% из статьи Grafke, Vanden-Eijnden, Nunmerical Computation of Rare Events
% via LDT. 2018. Chapter 6.C
% 
%  N - total population
% 
% 
% Measles - Корь, краснуха 
% Disease - болезнь, заболеваемость
% eradicate - искоренять, истреблять, ликвидировать
% 
% s - susceptible individuals (восприимчивая группа)
% i - infected individuals (инфицированные)
% r - recovery individuals, thus immune (выздоравлившие)
% v - vaccinated individuals, thus immune
% 
% individuals are born  and die with same rate Mu
% q - vaccination rate (probability of vaccinated or susceptible
% otherwise).
% Child is born vaccinated with probability q.
% 
%  Desease is transmitted by contact between i and s individuals with
%  contact rate Beta
%  Recovery rate is Gamma.
% 
% Ro is reproduction number, average number of contacts per infected
% individuals.
% Ro = beta/(mu+gamma)
% 
% If q >= q_thr, where q_thr = 1 - 1/Ro
% then disease eradicate eventually
% 
% The dynamic of s and i are independed of V and R.
% 


function [fright] = getVacModel2()
% особые точки
% ps = [gm/dl al/bt];

fright = @right;

end
%%

function dx = right(t,x, w, beta, mu, gamma, q, N)
    s = x(1);
    i = x(2);
    v = x(3);
    r = x(4);
    
    ds = mu*N*(1-q) - mu*s - beta*i*s/N + w(1);
    di = beta*s*i/N - (mu + gamma) * i + w(2);
    dv = mu*N*q - mu * v;
    dr = gamma*i + mu*r;
    
    dx = [ds; di; dv; dr];
end


