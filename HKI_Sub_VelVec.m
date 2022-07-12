function [Vu_Dimensionless]=HKI_Sub_VelVec(VelSurSel, u0, NN)
Nb = NN(1); Nr = NN(2); Nt = NN(3);

if VelSurSel == 1;
    % Bottom �鸸 1m/s�� ����
    Vu = zeros(Nb+Nr+Nt,1);
    Vu(1:Nb,1) = -1;
elseif VelSurSel == 2;
    % Cylinder �鸸 1m/s�� ����
    Vu = zeros(Nb+Nr+Nt,1);
    Vu(Nb+1:Nb+Nr,1) = 1;
elseif VelSurSel == 3;
    % Top �鸸 1m/s�� ����
    Vu = zeros(Nb+Nr+Nt,1);
    Vu(Nb+Nr+1:end,1) = 1;
elseif VelSurSel == 4;
    % ��� surface �� 1m/s �� ����
    Vu = ones(Nb+Nr+Nt,1);
    Vu(1:Nb,1) = -1;
end

%% Dimensionless
Vu_Dimensionless = Vu./u0;
end