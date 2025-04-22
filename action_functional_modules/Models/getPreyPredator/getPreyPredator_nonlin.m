%% ћодель хищник-жертва Ћотке и ¬ольтерра
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
% из ¬ики
% x=[x_1  x_2 ]^T Ц вектор состо€ни€, количество жертв и хищников соответственно, 
% ? Ц коэффициент рождаемости жертв, ? Ц коэффициент убыли хищников, 
% ? Ц коэффициент убийства жертв, ? Ц коэффициент воспроизводства сытых хищников, 
% g_1  и g_2 Ц параметры возмущений извне, ?(t) Ц белый шум (возмущение извне).

function [right, ps, right_migrt] = getPreyPredator_nonlin(al, bt, gm, dl, dtl)
% особые точки
ps = [gm/dl al/bt];

% нелинейна€ модель
right = @(t,x)[al*x(1) - bt*x(1)*x(2);    % количество жертв
           -gm*x(2) + dl*x(1)*x(2)] ;  % количество хищников

% from art. Numerical comp. p 11
% right = @(t,x)[al*x(1) - bt*x(1)*x(2) + dtl;    % количество жертв
%               -gm*x(2) + dl*x(1)*x(2) + dtl] ;  % количество хищников

       
right_migrt = @right_migration;

end
%%

% нелинейна€ модель с миграцей дл€ шума
% function dx = right_migration(t,x,w,wtime, al,bt,gm,dl,g1,g2)
% f = interp1(wtime,w,t, 'previous'); % Interpolate the data set (ft,f) at time t
% dx = zeros(2,1);
% dx(1) =  al*x(1) - bt*x(1)*x(2) + g1*f(1);    % количество жертв
% dx(2) = -gm*x(2) + dl*x(1)*x(2) + g2*f(2);         % количество хищников
% end

% function dx = right_migration(t,x,p)
% al = p.al;
% bt = p.bt;
% gm = p.gm;
% dl = p.dl;
% g1 = p.g1;
% g2 = p.g2;
% w = p.w;

% function dx = right_migration(t,x,  al,bt,gm,dl,g1,g2, w)
% dx = zeros(2,1);
% dx(1) =  al*x(1) - bt*x(1)*x(2) + g1*w(1);    % количество жертв
% dx(2) = -gm*x(2) + dl*x(1)*x(2) + g2*w(2);         % количество хищников
% end

% from art. Numerical comp. p 11
function dx = right_migration(t,x,  al,bt,gm,dl,g1,g2, w, dtl)
dx = zeros(2,1);
dx(1) =  al*x(1) - bt*x(1)*x(2) + dtl + g1*w(1);    % количество жертв
dx(2) = -gm*x(2) + dl*x(1)*x(2) + dtl + g2*w(2);         % количество хищников
end
