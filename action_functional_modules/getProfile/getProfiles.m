% Считаем несколько А-профилей
% 
% dT = 1/32;
% time = 0 : dT : 15;
% pfls = getProfiles(sys, time, thrs, eps);
% 
% 
% Изменения:
% 2022.02.17
%    Закомментировал вывод времени. Этот модуль понадоился в во внешнем
%    цикле, который много раз выполняется.


%%
function pfls = getProfiles(sys, time, thrs, eps)
pfls = [];

% fprintf('Считаем А-профили... '), tic
for thr = thrs
    pfl = getProfile_v2(sys, time, thr);
    pfl.pr = exp(-eps^-2 * pfl.aaf) ;
    pfl.prVf = exp(eps^-2 * pfl.vf) ;
    pfl.estTime = exp(eps^-2 * pfl.aaf);  % оценка среднего времени выхода к порогу
    pfls = [pfls pfl]; 
end
% toc

end