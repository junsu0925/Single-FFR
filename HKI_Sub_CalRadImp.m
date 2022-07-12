function [VF_For_FFR] = HKI_Sub_CalRadImp(InvMat_HKI,radius_a,length_l,NN)
Nb = NN(1); Nr = NN(2); Nt = NN(3);
TotNumPre = (Nb + Nr + Nt) - 2;

VA_b = zeros(1,Nb);
VA_r = zeros(1,Nr);

%% Area Vector define
VA_b(1,1) = 1/6;
VA_b(1,2:Nb-1) = (1:Nb-2);
VA_b(1,end) = (3*(Nb-1)-1)/6;
VA_b = -VA_b.*((2*pi)/(Nb-1)^2);

VA_t = -flipud(VA_b')';

VA_r(1,1) = 1/2;
VA_r(1,2:Nr-1) = 1;
VA_r(1,end) = 1/2;
VA_r = VA_r.*((2*pi)/(Nr-1));

%% Radiation Impedance Matrix Define (same as Force)

Nb = NN(1); Nr = NN(2); Nt = NN(3);
TotNumPre = (Nb + Nr + Nt) - 2;

T_3B_1 = zeros(Nb,3);
T_3B_1(:,1) = 1;
T_3B_2 = zeros(Nr,3);
T_3B_2(:,2) = 1;
T_3B_3 = zeros(Nt,3);
T_3B_3(:,3) = 1;
T_3B = [T_3B_1; T_3B_2; T_3B_3];

GA = zeros(3,TotNumPre);
VA_b = zeros(1,Nb);
VA_r = zeros(1,Nr);

%% Area Vector define
VA_b(1,1) = 1/6;
VA_b(1,2:Nb-1) = (1:Nb-2);
VA_b(1,end) = (3*(Nb-1)-1)/6;
VA_b = VA_b.*((2*pi)/(Nb-1)^2);

VA_t = flipud(VA_b')';

VA_r(1,1) = 1/2;
VA_r(1,2:Nr-1) = 1;
VA_r(1,end) = 1/2;
VA_r = VA_r.*((2*pi)/(Nr-1));

GA(1,1:Nb) = -VA_b; % Bottom Force, (-)
GA(2,Nb:Nb+Nr-1) = VA_r;
GA(3,Nb+Nr-1:end) = VA_t;

%% Radiation Impedance Matrix Define (same as Force)
HKI_to_FFR_AreaRatio_M = [radius_a/(2*pi*length_l) 0 0;0 1/(2*pi) 0;0 0 radius_a/(2*pi*length_l)];
VF_For_FFR = HKI_to_FFR_AreaRatio_M*GA*InvMat_HKI*T_3B;

end


