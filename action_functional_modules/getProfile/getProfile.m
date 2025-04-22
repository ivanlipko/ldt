%% ������ ��������� ������� � �������������� �������� ����������� ��������
%
% ���������:
% 19.08.2020
%     ������� �� ������ ������ �������
% 05.09.2019
%     ����������� pfl.af ������ pfl.aaf . ������ ��� ��� �-�������
% 15.06.2019
%     ������� ����� �������.
%     ��������� ������ ������� ��� ���������, ������� �������� �����,
%     ���������� � ����.
% 25.06.2019
%     �����������
% 26.10.2018
%     ����� �� ������������ ��������� �� ������, ������ ������ ������������.
% 19.10.2018
%     �����������
% 02.10.2018
%     ������� ��������� ��������� �� UTF-8
%     ������� ��������� ������� �� �������� ���������
% 17.04.2018
%     ������� ������������ �������� ������� ������� � � �������/������� �
%        ����������, ���� ��� ����������.
%     �����������
%
% % ��� ���������� ���������� - ����� ���� ��� ������� �����������

%% TODO
% - ����� ���� �����, ���� �� ���� ����� SS-������� sys, ����� � �.�.,
% ����������� ���������� ����� ������������, ����� ����� �����
% �������������� �������, � ������ ��������������

%%
% function [ssG_out, time, prfle, psiTf, ...
%     prfPeterb, covMatrix, sysA] = ...
%     getProfile(sys, sysG, xTf, kedTf, d_time, cols, ssPeterb)
% function [time, path, peterb, af] = ...
%     getProfile(sys, sysG, xTf, Tf, dT, cols, ssPeterb)
function [pfl] = ...
    getProfile(sys, sysG, xTf, Tf, dT, cols, ssPeterb)

ASize = size(sys.A,1);

%% ���������� �������������� �������
% ��� ����� ������ ������ ssModel � ������ ��� ��� ����� ������� ����������
% (�� ��������� ����� ���������)
% ��� ����� ��� ���������� ����� ����������, ������� �������� �������� ���
% ���������� ������� ������ ���������� � ����� xTf

if nargin == 7
    % ������ � ������ ������ ����������
    ssPeterbA = ssPeterb.A;
    ssPeterbC = ssPeterb.C;
    ssPeterbG = ssPeterb.G;
    sysAPeterbSize = size(ssPeterbA,1);
    sysA = [
        sys.A    sysG*ssPeterbC(1,:)
        zeros(sysAPeterbSize, ASize)   ssPeterbA  ];
    % ssG = [
    %     zeros(sizeofAssFeedbacked(1), 1)
    %     ssPeterbG];
    sysG = [ % ssG_out = ?
        sysG
        ssPeterbG];
    
    % ������������ ����� ������� ����������
    %  ���� ������ �������. ���� �� ���������������
    % ���������� ��� xTf ������ ���� � ������ ������. ����� ���� ����� ������
    %     ������������ ����� � ��������
    %         sysA = sysA([4,1,2,3,5,6,7,8,9,10,11],[4,1,2,3,5,6,7,8,9,10,11]);
    %         ssG_out = sysG([4,1,2,3,5,6,7,8,9,10,11]);
    ASize = ASize + sysAPeterbSize;  % �������� ������ ������� �
    
    % ��� �������� �������� ������
    %  ------------------------------------------------------------------------
    % [ sys.A    sysG*ssPeterbC(1,:) ];
    % [ zeros(sysAPeterbSize, sysASize-1)   ssPeterbA  ]    ;
    % [
    %     sys.A;%(1:4,1:4) ;
    %     zeros(sysAPeterbSize, sysASize )  ];%-1)  ];
    % [
    %     sysG_out*ssPeterbC(1,:)
    %     ssPeterbA  ];
    %  ------------------------------------------------------------------------
else
    % ������ ��� ����� ������ ����������
    
    if exist('cols', 'var')
        % ������������ ����� ������� ����������
        % ���������� ��� xTf ������ ���� � ������ ������.
        sysA = sys.A(cols,cols);
        sysG = sysG(cols, :);
    else
        %  ��� ������������, 17.04.2018
        sysA = sys.A;
    end
    % State of fullSystemA is
    % [
    %   ...
    % ]
end

covMat = lyap(sysA, sysG*sysG');

%% ���������� ������� ��� ���f
Dz = covMat(1,2:ASize);
D1 = covMat(1,1);
xTf = [xTf;  Dz'*xTf/D1];

% psiTf = inv(covMatrix) * xTf;
psiTf = covMat\ xTf;
% SF = psiTf' * covMatrix * psiTf; - �� ����� ���� ������

%% �������������. ���������� �������
% ������������� ����� ����� � �� ��������� ��������������� ��������� - ��� ���� ���� ��� �����������
% kedTime = kedTf - TT : 0.5*keddT : kedTf; % ������ -��..kedTf % ������������ ������. �� ���� ��� �� �����������. ������� ������� ����.
time = 0 : dT : Tf;
apath = zeros(ASize, length(time)); % SizeofFullSystemA-1
aaf = zeros( 1, length(time));
for i=1:length(time)
%     pfl(i,:) = expm( sysA * (time(i)-Tf) ) * (xTf - covMat*psiTf) + ...
%         + covMat * expm( sysA' * (Tf-time(i)) ) * psiTf;
    apath(:,i) = covMat * expm( sysA' * (Tf-time(i)) ) * psiTf;
    
    J0 =  covMat*expm((Tf-time(i))*sysA') - expm((time(i)-Tf)*sysA)*covMat; %(45)
    aaf(i) = -.5* psiTf'*expm((Tf-time(i))*sysA)*J0*psiTf; %(57);
end

ptb = zeros(size(sysG,2), length(time)); % SizeofFullSystemA-1
for i=1:length(time)
    ptb(:,i) =  sysG' * expm( sysA' * (Tf-time(i)) ) * psiTf;
end

%%
pfl.time = time;
pfl.apath = apath;
pfl.ptb = ptb;
pfl.aaf = aaf;

