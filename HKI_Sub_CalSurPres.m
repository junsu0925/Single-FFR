function [VP_nonDim, VP_Dim, Mat_3] = HKI_Sub_CalSurPres(GD, GW, p0, u0, NN)

%% Define Velocity Vector
IdentityGW = eye(length(GW));
Mat_1 = (IdentityGW - GW);
% Mat_2 = inv(Mat_1);
Mat_3 = Mat_1\GD;

% Surface 1, Bottom
[Vu]=HKI_Sub_VelVec(1, u0, NN);
Mat_4 = Mat_3*Vu;

VP_nonDim(:,1) = Mat_4;
VP_Dim(:,1) = VP_nonDim(:,1).*p0;

% Surface 2, Side
[Vu]=HKI_Sub_VelVec(2, u0, NN);
Mat_4 = Mat_3*Vu;

VP_nonDim(:,2) = Mat_4;
VP_Dim(:,2) = VP_nonDim(:,2).*p0;

% Surface 3, Top
[Vu]=HKI_Sub_VelVec(3, u0, NN);
Mat_4 = Mat_3*Vu;

VP_nonDim(:,3) = Mat_4;
VP_Dim(:,3) = VP_nonDim(:,3).*p0;

end