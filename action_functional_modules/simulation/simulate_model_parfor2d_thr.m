%% ������������� ������� �����-����� � ��������� ����������, ������������ �����
% ���������� ������ ������ ��� �����������

% ������ �������������: 
% -------------------------------------------------
% arWindSpeed = 10;
% sys = getKatWindModel(arWindSpeed, 0);
% disp(sys.StateName)
% dT = 0.1;
% 
% rng(123)
% exceedCount = 500;
% tFinal = 600;
% dT = 0.01;
% simTime = 0 : dT : tFinal;
% 
% xl0 = zeros(1, size(sys.A, 1));
% thr = 4
% 
% tic, fprintf('Simulation... ')
% exceed_outs = simulate_model_parfor2d_thr(sys, simTime, xl0, thr, exceedCount);
% exceedCount = size(exceed_outs,3)  % update count
% -------------------------------------------------
% 
% Plot simulation results:
% -------------------------------------------------
% figure(1), clf
% hold on, grid on
% xlabel(sys.StateName(1)), ylabel(sys.StateName(2)), title('������������ ������� ���������')
% 
% x = squeeze(exceed_outs(:,1,:));
% y = squeeze(exceed_outs(:,2,:));
% [n,c] = hist3([x(:) y(:)], 'Nbins',[200 200]);   % ��� ����� ���� ����� ����� ������
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
% ������������� � ������� ������ ����
% 
% https://www.gaussianwaves.com/2013/11/simulation-and-analysis-of-white-noise-in-matlab/
%
%%

function [exceed_outs,k] = simulate_model_parfor2d_thr(sys, time, x0, thr, exceedCount, batchSize)
% exceed_outs = zeros(length(time), 2, exceedCount);
exceed_outs = [];
% batchSize = exceedCount*2;
j = 0;
k = 0;
while j < exceedCount
    batch_outs = zeros(length(time), 2, batchSize);
%     batch_inds = zeros(batch, 1);
    batch_inds = [];
    
    parfor i=1:batchSize
        w = randn(length(time), size(sys.B,2));
        out = lsim(sys, w, time, x0);
        batch_outs(:,:,i) = out(:,1:2);

        % ����� ������ �����
        if thr > 0
            exceeds_mask = squeeze( out(:,1) >= thr );
        else
            exceeds_mask = squeeze( out(:,1) <= thr );
        end
        
        % �������� ���, ��� ����� �� �����
        if sum(exceeds_mask) > 0
            batch_inds = [batch_inds i];
        end
    end

    exceed_outs = cat(3,exceed_outs, batch_outs(:,:,batch_inds));

    j = j + length(batch_inds);
    fprintf('%d / %d \n', j, exceedCount);
    k = k + batchSize;
end

% fprintf('����� ������������� %d \n', k);

end

