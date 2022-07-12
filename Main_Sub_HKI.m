function [InvMat_HKI,VF_For_FFR,p0,Time_Main_Sub_HKI] = Main_Sub_HKI(totfreq,u0,c0,rho0,radius_a,length_l,NN)
%%
tic

%%
VF_nonDim = zeros(length(totfreq),5); VF_Dim = zeros(length(totfreq),5);
InvMat_HKI = cell(length(totfreq),1); VF_For_FFR = cell(length(totfreq),1);

handler = waitbar(0,'Initializing waitbar...'); % waitbar�� ���ϴ�.
for NumFreq = 1:length(totfreq)
    freq = totfreq(NumFreq);
    waitbar(NumFreq/length(totfreq),handler,sprintf('Computing... %d Hz', freq));
        
     [GD, GW, p0] = HKI_Sub_SurPres_Mat(freq,u0,c0,rho0,radius_a,length_l,NN);
    
     [VP_nonDim, VP_Dim, InvMat_HKI{NumFreq,1}] = HKI_Sub_CalSurPres(GD, GW, p0, u0, NN);
        
     [VF_For_FFR{NumFreq,1}] = HKI_Sub_CalRadImp(InvMat_HKI{NumFreq,1},radius_a,length_l,NN);

end
close(handler) % ������ ������ waitbar�� �����Ѵ�.

Time_Main_Sub_HKI = toc;
Hour = fix(Time_Main_Sub_HKI/3600); Min = fix(rem(Time_Main_Sub_HKI,3600)/60); Sec = round(rem(rem(Time_Main_Sub_HKI,3600),60));
fprintf('Main_Sub_HKI ��� �ҿ� �ð��� %d�ð� %d�� %d�� �Դϴ�.\n',Hour,Min,Sec)
end