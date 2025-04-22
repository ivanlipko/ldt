% �������� ���������� �� ������
% ���������� ���������� ����  �� �������� �� ������. � �������, ��� ���� ������ ������� ������� � ����� ������ �������� �� ������
% � ������ �� ��������� � ����� ������ ����� � ������ ���������, �� ���� ������������� �����������

% inSignal = Time X State X SimCount

% function [outs] = simToThr(nsys, xstate0, time, inSignal, right, simCount)
function [exceeds] = simToThr(thr, nsys, xstate0, time, peterb, rightfun, allCount, batchCount, time_thr)
dt = time(2) - time(1);
exceeds.paths = [];
exceeds.ptb = [];
exceeds.freq = 0;

tic, fprintf('������... ')
k = 0;   % ������� ���������
j = 0;   % ������� ��������� ������ ����������
l = 0;  % ������� ������
while j < allCount
    batch_outs = NaN(length(time), length(xstate0), batchCount);
    batch_w = NaN(length(time), size(peterb,2), batchCount);
    %     batch_inds = zeros(batch, 1);
    batch_inds = [];
    
    parfor c=1:batchCount
        state = xstate0;
        out = NaN(length(time), nsys);  %  ������������� �������
        z = 0;   % �������, ���������� �������. ��� ��� ������������ ��������
        for i = 1:length(time)
            out(i,:) = state;    % ���������� �.�.
            state = rkStep(time(i), dt, state, rightfun, peterb(i,:,c));
            if state(1) > thr    % ���� �������� ������
                if time(i)>time_thr   % ����� ��������� ������?
                    out(i,:) = state;   % ��������� ��������� ���������
                    out = circshift(out, length(time)-i);  % ���������� ������, ����� ����� ��� � ������ ����
                    batch_outs(:,:,c) = out;
                    batch_w(:,:,c) = peterb(:,:,c);
                    batch_inds = [batch_inds c];
                    break    % �����������
                end
                z = z+1;
            end
        end
    end    

    exceeds.paths = cat(3, exceeds.paths, batch_outs(:,:,batch_inds));
    exceeds.ptb = cat(3, exceeds.ptb, batch_w(:,:,batch_inds));    
   
    j = j + length(batch_inds);
    fprintf('%d / %d \n', j, allCount);
    k = k + batchCount;
    
    % ���������� ����������, ������ ��� �� ���������� ��������� �����
    if l > 10  
        break
    end
    l = l+1;
end
exceeds.freq = j/k;
toc
end



%% ���� ��� ������ �����-����� 4-�� �������
function f = rkStep(t, h, X, right, inSignal)

h2 = 0.5*h;
h6 = 0.166666666*h;
Fs = right(t,X, inSignal);
t = t + h2;

Xr = X + h2 * Fs;

F = right(t,Xr, inSignal);

s=F;
Fs=Fs+s+s;
Xr=X+h2*s;

F = right(t,Xr, inSignal);
t=t+h2;

s=F;
Fs=Fs+s+s;
Xr=X+h*s;

F = right(t,Xr, inSignal);

X=X+h6*(Fs+F);

f=X;
end
