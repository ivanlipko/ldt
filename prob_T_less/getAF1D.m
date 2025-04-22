% Вычисление ФД по первой данной от ВФ формуле 
% S = int(X)[ _t0^tf ( abs(dotX - b(X))^2  ) ds]

function af = getAF1D(sys, path)
dx = [ zeros(1,size(path,2)); diff(path)];
% dx = [ diff(path); zeros(1,size(path,2))];
% af = 0;
% 
% for i = 1:size(path,1)
%     x = path(i);
%     d = dx(i);
%     
%     b = sys.A*x;  % b(x)
%     q = (d - b)' * (d - b);
%     af = af + q;
% end
% 
% af = 0.5 * af;


bb = sys.A * path;
af = (dx - bb)' * (dx - bb)/2;

end
