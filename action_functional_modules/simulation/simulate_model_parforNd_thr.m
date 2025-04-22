%% Симулирование методом Монте-Карло и отcеивание траекторий, пересекающих порог
% возвращаем только заданные состояния statevec
% 
% Пример использования:
% -------------------------------------------------
% arWindSpeed = 10;
% sys = getKatWindModel(arWindSpeed, 0);
% disp(sys.StateName)
% dT = 0.1;
%
% rng(123)
% mmcCount = 1000; batchSize = 5000;
% tFinal = 600;
% dT = 0.01;
% simTime = 0 : dT : tFinal;
%
% x0 = zeros(1, size(sys.A, 1));
% thr = 4
% statevec = 1 : size(sys.A, 1);
%
% tic, fprintf('Simulation... ')
% [mmc_outs, k, mmc_w] = simulate_model_parforNd_thr(sys, simTime, x0, thr, statevec, mmcCount, batchSize);
% mmcCount = size(mmc_outs,3);  % update count
% toc
% -------------------------------------------------
%
% Plot simulation results:
% -------------------------------------------------
% figure(1), clf
% hold on, grid on
% xlabel(sys.StateName(1)), ylabel(sys.StateName(2)), title('Посещаемость фазовой плоскости')
%
% x = squeeze(exceed_outs(:,1,:));
% y = squeeze(exceed_outs(:,2,:));
% [n,c] = hist3([x(:) y(:)], 'Nbins',[200 200]);   % вот здесь надо очень много памяти
% n = n';
% cx = c{1};
% cy = c{2};
%
% pcolor(cx,cy, log(n+1))
% shading flat
% c = colorbar;
% c.Label.String = 'log(n+1)';

% -------------------------------------------------
%
% save 'simulation.mat' outs w time
%
%
% Симулирование с помощью белого шума
%
% https://www.gaussianwaves.com/2013/11/simulation-and-analysis-of-white-noise-in-matlab/
%
%%
% 2021.07.04 добавлена возможность использовать свою функцию для
% симулирования sfun
% 

function [exceed_outs, k, exceed_w] = simulate_model_parforNd_thr(...
    sys, time, x0, thr, statevec, mmcCount, batchSize, sfun)

if nargin == 7
    sfun = @simfun;
else
end

exceed_outs = [];
exceed_w = [] ;

j = 0;
k = 0;
while j < mmcCount
    batch_outs = zeros(length(time), length(statevec), batchSize);
    batch_w = zeros(length(time), size(sys.B,2), batchSize);
    %     batch_inds = zeros(batch, 1);
    batch_inds = [];
    
    
    parfor i=1:batchSize
%     for i=1:batchSize
        [out, w] = sfun(sys, time, x0);

        batch_outs(:,:,i) = out(:,statevec);
        batch_w(:,:,i) = w;
        
        % маска выбора точек
        if thr > 0
            exceeds_mask = squeeze( out(:,1) >= thr );
        else
            exceeds_mask = squeeze( out(:,1) <= thr );
        end
        
        % помечаем тех, кто вышел за порог
        if sum(exceeds_mask) > 0
            batch_inds = [batch_inds i];
        end
    end
    
    exceed_outs = cat(3, exceed_outs, batch_outs(:,:,batch_inds));
    exceed_w = cat(3, exceed_w, batch_w(:,:,batch_inds));
    
    j = j + length(batch_inds);
    fprintf('%d / %d \n', j, mmcCount);
    k = k + batchSize;
end

fprintf('Всего смоделировано %d \n', k);

end


% Функция симулирования 
function [out, w] = simfun(sys, time, x0)
w = randn(length(time), size(sys.B,2));
out = lsim(sys, w, time, x0);
end