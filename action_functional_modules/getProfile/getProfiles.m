% ������� ��������� �-��������
% 
% dT = 1/32;
% time = 0 : dT : 15;
% pfls = getProfiles(sys, time, thrs, eps);
% 
% 
% ���������:
% 2022.02.17
%    ��������������� ����� �������. ���� ������ ���������� � �� �������
%    �����, ������� ����� ��� �����������.


%%
function pfls = getProfiles(sys, time, thrs, eps)
pfls = [];

% fprintf('������� �-�������... '), tic
for thr = thrs
    pfl = getProfile_v2(sys, time, thr);
    pfl.pr = exp(-eps^-2 * pfl.aaf) ;
    pfl.prVf = exp(eps^-2 * pfl.vf) ;
    pfl.estTime = exp(eps^-2 * pfl.aaf);  % ������ �������� ������� ������ � ������
    pfls = [pfls pfl]; 
end
% toc

end