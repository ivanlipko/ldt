%% Модель хищник-жертва Лотке и Вольтерра
% линейная система (https://ru.wikipedia.org/wiki/Модель Лотки — Вольтерры)
% x0 = [ps(1) ps(2)] - x0;
% 
% al = 4;
% bt = .05;
% gm = 2;
% dl = .05;
% 
% al = 1;
% bt = .14;
% gm = .42;
% dl = .07;
%
%% Экспорт в word:
% https://www.youtube.com/watch?v=3jfLmmcG8fw     на  10 минуте
%             \matrix(0 & -bt*gm/dl @ dl*al/bt & 0)
% нажать пробел здесь^               а затем здесь ^
% 

%%
function [sys, ps] = getPreyPredator_linear(al, bt, gm, dl, g1, g2)
% особые точки
ps = [gm/dl al/bt];

A = [0 -bt*gm/dl;
     dl*al/bt 0];

G = [ g1  0     
      0   g2];  
C = eye(size(A));
D = zeros(size(G));

sys = ss(A,G,C,D);
sys.Name = 'Хищник-Жертва (Лотке-Вольтерра). Модель в отклонениях от стационарной точки';
sys.StateName = {'Жертвы', 'Хищники' };
sys.StateUnit = {'шт','шт'};
sys.InputName = {'Жертвы-миграция', 'Хищники-миграция'};
sys.OutputName = sys.StateName;
end

