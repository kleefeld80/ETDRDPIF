% Experiment to compute rate of convergence of ETD RRP in comparison with
% other second ETD schemes
% November 2014

clc;clear
%k1= 0.1;k2=0.05;k3=0.025;k4=0.0125;k5=0.00625;
k1=0.05;k2=0.025;k3=0.0125;k4=0.00625;k5=0.003125;
%k1=0.025;k2=0.0125;k3=0.00625;k4=0.003125;k5=0.001563;

%h1 = 9; h2 = 19; h3=39; h4=79; h5=159;
h1=79;h2=h1;h3=h1;h4=h1;h5=h1;
step = [k1,k2,k3,k4];
space=[h1,h2,h3,h4];
%save step 
%% results for BDF2
[runtime1,soln1] = enzymekinetics2Dbdf2(k1,h1,k1^3);
[runtime2,soln2] = enzymekinetics2Dbdf2(k2,h2,k2^3);
[runtime3,soln3] = enzymekinetics2Dbdf2(k3,h3,k3^3);
[runtime4,soln4] = enzymekinetics2Dbdf2(k4,h4,k4^3);
[~,soln5] = enzymekinetics2Dbdf2(k5,h5,k5^3);

 solnref2 = soln2; solnref3 = soln3; solnref4 = soln4; solnref5 = soln5;
 TimeBDF2=[runtime1,runtime2,runtime3,runtime4];
 [convBDF2,errorBDF2]=conv_ekinetics(soln1,soln2,soln3,soln4,solnref2,solnref3,solnref4,solnref5,step);
 
 save EkineticsDATA0202_smooth_79_fine solnref2 solnref3 solnref4 solnref5 TimeBDF2 convBDF2 errorBDF2
 %load EkineticsDATA0202

%% results for RDP 
[runtime1,soln1] = enzymekinetics_2D_ETDRDP(k1,h1);
[runtime2,soln2] = enzymekinetics_2D_ETDRDP(k2,h2);
[runtime3,soln3] = enzymekinetics_2D_ETDRDP(k3,h3);
[runtime4,soln4] = enzymekinetics_2D_ETDRDP(k4,h4);
%[runtime5,soln5] = enzymekinetics_2D_ETDRDP(k5,h5);

 %solnref2 = soln2; solnref3 = soln3; solnref4 = soln4; solnref5 = soln5;
 TimeRDP=[runtime1,runtime2,runtime3,runtime4];
 [convRDP,errorRDP]=conv_ekinetics(soln1,soln2,soln3,soln4,solnref2,solnref3,solnref4,solnref5,step);

%% results for RDP-split 
[runtime1,soln1] = enzymekinetics_2D_ETDRDP_split(k1,h1);
[runtime2,soln2] = enzymekinetics_2D_ETDRDP_split(k2,h2);
[runtime3,soln3] = enzymekinetics_2D_ETDRDP_split(k3,h3);
[runtime4,soln4] = enzymekinetics_2D_ETDRDP_split(k4,h4);
%[runtime5,soln5] = enzymekinetics_2D_ETDRDP_split(k5,h5);

 %solnref2 = soln2; solnref3 = soln3; solnref4 = soln4; solnref5 = soln5;
 TimeRDPS=[runtime1,runtime2,runtime3,runtime4];
 [convRDPS,errorRDPS]=conv_ekinetics(soln1,soln2,soln3,soln4,solnref2,solnref3,solnref4,solnref5,step); 
 
% %% results for ETDCN
% [runtime1,soln1] = enzymekinetics_2D_ETDCN(k1,h1);
% [runtime2,soln2] = enzymekinetics_2D_ETDCN(k2,h2);
% [runtime3,soln3] = enzymekinetics_2D_ETDCN(k3,h3);
% [runtime4,soln4] = enzymekinetics_2D_ETDCN(k4,h4);
% %[~,soln5] = enzymekinetics_2D_ETDCN(k5,h5);
% 
%  TimeCN=[runtime1,runtime2,runtime3,runtime4];
%  [convCN,errorCN]=conv_ekinetics(soln1,soln2,soln3,soln4,solnref2,solnref3,solnref4,solnref5,step);
% 
%  %% results for Pade02
% [runtime1,soln1] = enzymekinetics_2D_ETDP02(k1,h1);
% [runtime2,soln2] = enzymekinetics_2D_ETDP02(k2,h2);
% [runtime3,soln3] = enzymekinetics_2D_ETDP02(k3,h3);
% [runtime4,soln4] = enzymekinetics_2D_ETDP02(k4,h4);
% %[~,solnref] = enzymekinetics_2D_ETDP02(k5,h5);
% 
%  TimeP02=[runtime1,runtime2,runtime3,runtime4];
%  [convP02,errorP02]=conv_ekinetics(soln1,soln2,soln3,soln4,solnref2,solnref3,solnref4,solnref5,step);
 %% results for IMEX-BDF2
[runtime1,soln1] = enzymekinetics_2D_IMEXBDF2(k1,h1);
[runtime2,soln2] = enzymekinetics_2D_IMEXBDF2(k2,h2);
[runtime3,soln3] = enzymekinetics_2D_IMEXBDF2(k3,h3);
[runtime4,soln4] = enzymekinetics_2D_IMEXBDF2(k4,h4);
%[runtime5,soln5] = enzymekinetics_2D_IMEXBDF2(k5,h5);

% [runtime1,soln1] = enzymekinetics_2D_IMEXBDF2_v2(k1,h1);
% [runtime2,soln2] = enzymekinetics_2D_IMEXBDF2_v2(k2,h2);
% [runtime3,soln3] = enzymekinetics_2D_IMEXBDF2_v2(k3,h3);
% [runtime4,soln4] = enzymekinetics_2D_IMEXBDF2_v2(k4,h4);
% [runtime5,soln5] = enzymekinetics_2D_IMEXBDF2_v2(k5,h5);

 %solnref2 = soln2; solnref3 = soln3; solnref4 = soln4; solnref5 = soln5;

 TimeIM2=[runtime1,runtime2,runtime3,runtime4];
 [convIM2,errorIM2]=conv_ekinetics(soln1,soln2,soln3,soln4,solnref2,solnref3,solnref4,solnref5,step);
  %% results for IMEX-TR
[runtime1,soln1] = enzymekinetics_2D_IMEXTR(k1,h1);
[runtime2,soln2] = enzymekinetics_2D_IMEXTR(k2,h2);
[runtime3,soln3] = enzymekinetics_2D_IMEXTR(k3,h3);
[runtime4,soln4] = enzymekinetics_2D_IMEXTR(k4,h4);
%[runtime5,soln5] = enzymekinetics_2D_IMEXTR(k5,h5);
 %solnref2 = soln2; solnref3 = soln3; solnref4 = soln4; solnref5 = soln5;

 TimeIMTR=[runtime1,runtime2,runtime3,runtime4];
 [convIMTR,errorIMTR]=conv_ekinetics(soln1,soln2,soln3,soln4,solnref2,solnref3,solnref4,solnref5,step);
  %% results for IMEX-CNLF
[runtime1,soln1] = enzymekinetics_2D_IMEXCNLF(k1,h1);
[runtime2,soln2] = enzymekinetics_2D_IMEXCNLF(k2,h2);
[runtime3,soln3] = enzymekinetics_2D_IMEXCNLF(k3,h3);
[runtime4,soln4] = enzymekinetics_2D_IMEXCNLF(k4,h4);
%[runtime5,soln5] = enzymekinetics_2D_IMEXCNLF(k5,h5);
 %solnref2 = soln2; solnref3 = soln3; solnref4 = soln4; solnref5 = soln5;

 TimeIMCNLF=[runtime1,runtime2,runtime3,runtime4];
 [convIMCNLF,errorIMCNLF]=conv_ekinetics(soln1,soln2,soln3,soln4,solnref2,solnref3,solnref4,solnref5,step);
 
  %% results for IMEX-Adams2
[runtime1,soln1] = enzymekinetics_2D_IMEXAdams2(k1,h1);
[runtime2,soln2] = enzymekinetics_2D_IMEXAdams2(k2,h2);
[runtime3,soln3] = enzymekinetics_2D_IMEXAdams2(k3,h3);
[runtime4,soln4] = enzymekinetics_2D_IMEXAdams2(k4,h4);
%[runtime5,soln5] = enzymekinetics_2D_IMEXAdams2(k5,h5);
 %solnref2 = soln2; solnref3 = soln3; solnref4 = soln4; solnref5 = soln5;

 TimeIMAD2=[runtime1,runtime2,runtime3,runtime4];
 [convIMAD2,errorIMAD2]=conv_ekinetics(soln1,soln2,soln3,soln4,solnref2,solnref3,solnref4,solnref5,step);
%% Produce Efficiency plot
% ETD Comparison
% Time_mat = [TimeRDP;TimeCN;TimeP02;TimeIM2];
% Error_mat = [errorRDP;errorCN;errorP02;errorIM2];

% IMEX Comparison
Time_mat = [TimeRDP;TimeIM2;TimeIMTR;TimeIMCNLF;TimeIMAD2;TimeRDPS];
Error_mat = [errorRDP;errorIM2;errorIMTR;errorIMCNLF;errorIMAD2;errorRDPS];
% Time_mat = [TimeRDP;TimeIM2;TimeIMTR;TimeIMCNLF;TimeIMAD2];
% Error_mat = [errorRDP;errorIM2;errorIMTR;errorIMCNLF;errorIMAD2];


save EkineticsDATA0202 Time_mat Error_mat

efficiency_plot_ekinetics(Time_mat,Error_mat)
convergence_plot_ekinetics(step,Error_mat)
%% Display results
fprintf(' Results for ETDRDP\n')
fprintf('k             h            error        conv       Time \n');
for i = 1:4
    fprintf('%.6f   %.6f    %1.4e      %.2f      %.5f\n', step(i),1/(space(i)+1),errorRDP(i),convRDP(i),TimeRDP(i))
end
fprintf(' Results for ETDRDPsplit\n')
fprintf('k             h            error        conv       Time \n');
for i = 1:4
    fprintf('%.6f   %.6f    %1.4e      %.2f      %.5f\n', step(i),1/(space(i)+1),errorRDPS(i),convRDPS(i),TimeRDPS(i))
end
% fprintf('\n Results for ETDCN\n')
% fprintf('k             h            error        conv       Time \n');
% for i = 1:4
%     fprintf('%.6f   %.6f    %1.4e      %.2f      %.5f\n', step(i),1/(space(i)+1),errorCN(i),convCN(i),TimeCN(i))
% end
% 
% fprintf('\n Results for ETDP02\n')
% fprintf('k             h            error        conv       Time \n');
% for i = 1:4
%     fprintf('%.6f   %.6f    %1.4e      %.2f      %.5f\n', step(i),1/(space(i)+1),errorP02(i),convP02(i),TimeP02(i))
% end
fprintf('\n Results for IMEX-BDF2\n')
fprintf('k             h            error        conv       Time \n');
for i = 1:4
    fprintf('%.6f   %.6f    %1.4e      %.2f      %.5f\n', step(i),1/(space(i)+1),errorIM2(i),convIM2(i),TimeIM2(i))
end
fprintf('\n Results for IMEX-TR\n')
fprintf('k             h            error        conv       Time \n');
for i = 1:4
    fprintf('%.6f   %.6f    %1.4e      %.2f      %.5f\n', step(i),1/(space(i)+1),errorIMTR(i),convIMTR(i),TimeIMTR(i))
end
fprintf('\n Results for IMEX-CNLF\n')
fprintf('k             h            error        conv       Time \n');
for i = 1:4
    fprintf('%.6f   %.6f    %1.4e      %.2f      %.5f\n', step(i),1/(space(i)+1),errorIMCNLF(i),convIMCNLF(i),TimeIMCNLF(i))
end
fprintf('\n Results for IMEX-Adams2\n')
fprintf('k             h            error        conv       Time \n');
for i = 1:4
    fprintf('%.6f   %.6f    %1.4e      %.2f      %.5f\n', step(i),1/(space(i)+1),errorIMAD2(i),convIMAD2(i),TimeIMAD2(i))
end
