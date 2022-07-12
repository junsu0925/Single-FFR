clear all
close all

%% HKI Data Load
load('matlab.mat')

%% COMSOL Data Load
opts = delimitedTextImportOptions("NumVariables", 3);
opts.DataLines = [9, Inf];
opts.Delimiter = " ";
opts.VariableNames = ["VarName1", "Model", "Var3"];
opts.SelectedVariableNames = ["VarName1", "Model"];
opts.VariableTypes = ["double", "double", "string"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts.ConsecutiveDelimitersRule = "join";
opts.LeadingDelimitersRule = "ignore";
opts = setvaropts(opts, "Var3", "WhitespaceRule", "preserve");
opts = setvaropts(opts, "Var3", "EmptyFieldRule", "auto");
tbl = readtable("D:\OneDrive - postech.ac.kr\Graduate\X - 연구\2. 2D FFR\Single FFR Debug\Single FFR with HKI\COMSOL Data 비교\Conductance - Single FFR.txt", opts);
Freq = tbl.VarName1;
ConductanceSingleFFR = tbl.Model;
clear opts

opts = delimitedTextImportOptions("NumVariables", 2);
opts.DataLines = [9, Inf];
opts.Delimiter = " ";
opts.VariableNames = ["VarName1", "Model"];
opts.VariableTypes = ["double", "double"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts.ConsecutiveDelimitersRule = "join";
opts.LeadingDelimitersRule = "ignore";
tbl = readtable("D:\OneDrive - postech.ac.kr\Graduate\X - 연구\2. 2D FFR\Single FFR Debug\Single FFR with HKI\COMSOL Data 비교\TVR - Single FFR.txt", opts);
Freq = tbl.VarName1;
TVRSingleFFR = tbl.Model;
clear opts

f0 = 2.806030228019852e+03;
%% Plot
figure(1)
plot(ExportCircuit(:,2),real(ExportCircuit(:,3)),'k--','LineWidth',2)
hold on
plot(Freq/f0,ConductanceSingleFFR,'k','LineWidth',2)
hold on
grid on
xlabel('ka/k_0 a','fontsize',20, 'fontangle','italic');
ylabel('Conductance','fontsize',20, 'fontangle','italic');
set(gca, 'fontsize',16)
set(gcf, 'color', 'w')
legend('FFR - HKI','FFR - COMSOL')
xlim([min(Freq/f0) max(Freq/f0)])

figure(2)
plot(ExportCircuit(:,2),ExportObPres(:,3),'k--','LineWidth',2)
hold on
plot(Freq/f0,TVRSingleFFR,'k','LineWidth',2)
hold on
grid on
xlabel('ka/k_0 a','fontsize',20, 'fontangle','italic');
ylabel('TVR [dB]','fontsize',20, 'fontangle','italic');
set(gca, 'fontsize',16)
set(gcf, 'color', 'w')
legend('FFR - HKI','FFR - COMSOL')
xlim([min(Freq/f0) max(Freq/f0)])