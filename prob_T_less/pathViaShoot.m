%% Поиск профиля между двумя заданными точки

function [shoots] = pathViaShoot(time, sys, nu, y0, yf)
% fprintf('pathViaShoot ... ')

% Первый выстрел
pf0 = getCondition(nu(1));
[control] = calcPsolution(sys, pf0, time);
[sysOut, e1] = errSys(sys, control, time, y0, yf);

shoots.errors = [e1];
shoots.control = [control(:,1)];
shoots.yf = [sysOut(end,1)];
shoots.nu = [nu(1)];
shoots.syspath(1,:,:) = sysOut;

% второй выстрел
pf0 = getCondition(nu(2));
[control] = calcPsolution(sys, pf0, time);
[sysOut, e2] = errSys(sys, control, time, y0, yf);

shoots.errors = [shoots.errors; e2];
shoots.control = [shoots.control control(:,1)];
shoots.yf = [shoots.yf; sysOut(end,1)];
shoots.nu = [shoots.nu; nu(2)];
shoots.syspath(2,:,:) = sysOut;

eps = 10^-5;
max_iter = 7;
iter = 1;

% while(abs(e2-e1) > eps)
while(abs(nu(2)-nu(1)) > eps)
    nu2_old = nu(2);
    nu(2) = nu(2) - (e2 * (nu(2)-nu(1)))/(e2-e1);
    nu(1) = nu2_old;
    e1 = e2;
    
    pf0 = getCondition(nu(2));
    [control] = calcPsolution(sys, pf0, time);
    [sysOut, e2] = errSys(sys, control, time, y0, yf);
    
    shoots.errors = [shoots.errors; e2];
    shoots.control = [shoots.control control(:,1)];
    shoots.yf = [shoots.yf; sysOut(end,1)];
    shoots.nu = [shoots.nu; nu(2)];
    shoots.syspath(iter+2,:,:) = sysOut;
    
    if(iter>max_iter)
        break
    end
    iter = iter + 1;
    
%     fprintf('.%d',iter)
end
% fprintf('Finished\n')
% fprintf('%d ..Finished\n',iter)

end

function [yf] = getCondition(nu)
yf = [nu];
end


function [out, e] = errSys(sys, control, time, y0, yf)
[out,~] = lsim(sys, control, time, y0);
e = out(end,1) - yf(1);
end


% Решение ОДУ
function [control] = calcPsolution(sys, pf, time)
% [~,p] = ode45(@rSysP, fliplr(time), pf); % справа налево
% [~,p] = ode45(@rSysP, time, pf);

A = sys.A;
B = sys.B;

p = zeros(length(time),size(B,1));
for i=1:length(time)
    %     p(i,:) = expm( -A' * (time(i)) ) * pf ;
    p(i,:) = expm( -A' * (time(length(time)-i+1)) ) * pf ;   % справа налево
end

control = B' * p';
control = fliplr(control)';

% figure(100), plot(time, p)
end


function dx = rSysP(t, x, sys)
% dp = - A' * dx .*nu' + A' * A * x .*nu' - A' * p;  % критерий S = \int norm2(dot X - f(X))
% dp = - A' * dx + A' * A * x - A' * p;
dx = -sys.A' * x;  % критерий S = \int (v^T * v)
end

