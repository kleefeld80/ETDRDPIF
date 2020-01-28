% Experiment to compare rate of convergence of ETDRRP and ETDPADE02 and
% ETDCN for Brusselator in 2D
% November 2014


clc;clear
% Initial parameters
 k1= 0.1;k2=0.05;k3=0.025;k4=0.0125;kref=0.00625;
  %k1=0.05;k2=0.025;k3=0.0125;k4=0.00625;kref=0.003125;
%k1=0.025;k2=0.0125;k3=0.00625;k4=0.003125; kref=0.001563;
%k1=0.0125;k2=0.00625;k3=0.003125; k4 = 0.001563; kref = 0.000782;


 tol = 0.01;
 %h1 = 11; h2 = 21; h3 = 41; h4= 81; href = 161;
 h1 = 81; h2 = h1; h3=h1; h4=h1;href = h1;
 step = [k1,k2,k3,k4];
 space=[1/(h1-1),1/(h2-1),1/(h3-1),1/(h4-1),1/(href-1)];
 
% 
%results for BDF2 
[runtime1,soln1] = brusselator2Dbdf2(k1,h1,k1^3);
[runtime2,soln2] = brusselator2Dbdf2(k2,h2,k2^3);
[runtime3,soln3] = brusselator2Dbdf2(k3,h3,k3^3);
[runtime4,soln4] = brusselator2Dbdf2(k4,h4,k4^3);
[~,soln5] = brusselator2Dbdf2(kref,href,kref^3);

 solnref2 = soln2; solnref3 = soln3; solnref4=soln4; solnref5=soln5;
 TimeBDF2=[runtime1,runtime2,runtime3,runtime4];
 [convBDF2,errorBDF2]=conv_Brusselator(soln1,soln2,soln3,soln4,solnref2,solnref3,solnref4,solnref5,step);
 
 
 save BrusselatorDATA_81 solnref2 solnref3 solnref4 solnref5 TimeBDF2 convBDF2 errorBDF2
 %load BrusselatorDATA
 
%% results for ETDRDP
[runtime1,soln1] = Brusselator2D_ETDRDP(k1,h1);
[runtime2,soln2] = Brusselator2D_ETDRDP(k2,h2);
[runtime3,soln3] = Brusselator2D_ETDRDP(k3,h3);
[runtime4,soln4] = Brusselator2D_ETDRDP(k4,h4);
%[~,soln5] = Brusselator2D_ETDRDP(kref,href);

%solnref2 = soln2; solnref3 = soln3; solnref4=soln4; solnref5=soln5;
 TimeRDP=[runtime1,runtime2,runtime3,runtime4];
 [convRDP,errorRDP]=conv_Brusselator(soln1,soln2,soln3,soln4,solnref2,solnref3,solnref4,solnref5,step);
 
%% results for ETDRDPsplit
[runtime1,soln1] = BrusselatorLOD2D_IFETDRDP(k1,h1);
[runtime2,soln2] = BrusselatorLOD2D_IFETDRDP(k2,h2);
[runtime3,soln3] = BrusselatorLOD2D_IFETDRDP(k3,h3);
[runtime4,soln4] = BrusselatorLOD2D_IFETDRDP(k4,h4);
%[~,soln5] = BrusselatorLOD2D_IFETDRDP(kref,href);

%solnref2 = soln2; solnref3 = soln3; solnref4=soln4; solnref5=soln5;
 TimeRDPS=[runtime1,runtime2,runtime3,runtime4];
 [convRDPS,errorRDPS]=conv_Brusselator(soln1,soln2,soln3,soln4,solnref2,solnref3,solnref4,solnref5,step); 
 
%% results for ETDCN
% [runtime1,soln1] =  Brusselator2D_ETDCN(k1,h1);
% [runtime2,soln2] =  Brusselator2D_ETDCN(k2,h2);
% [runtime3,soln3] = Brusselator2D_ETDCN(k3,h3);
% [runtime4,soln4] = Brusselator2D_ETDCN(k4,h4);
% [~,solnref] = Brusselator2D_ETDCN(kref,href);
% 
%  TimeCN=[runtime1,runtime2,runtime3,runtime4];
%  [convCN,errorCN]=conv_Brusselator(soln1,soln2,soln3,soln4,solnref2,solnref3,solnref4,solnref5,step);


%% results for ETDPADE02
% [runtime1,soln1] =  Brusselator2D_ETDpade02(k1,h1);
% [runtime2,soln2] =  Brusselator2D_ETDpade02(k2,h2);
% [runtime3,soln3] = Brusselator2D_ETDpade02(k3,h3);
% [runtime4,soln4] = Brusselator2D_ETDpade02(k4,h4);
% %[~,solnref] = Brusselator2D_ETDpade02(kref,href);
% 
%  Timepa=[runtime1,runtime2,runtime3,runtime4];
%  [convpa,errorpa]=conv_Brusselator(soln1,soln2,soln3,soln4,solnref2,solnref3,solnref4,solnref5,step);

%% results for IMEX-BDF2
[runtime1,soln1] =  Brusselator2D_IMEX_BDF2(k1,h1);
[runtime2,soln2] =  Brusselator2D_IMEX_BDF2(k2,h2);
[runtime3,soln3] = Brusselator2D_IMEX_BDF2(k3,h3);
[runtime4,soln4] = Brusselator2D_IMEX_BDF2(k4,h4);
%[~,solnref] = Brusselator2D_ETDpade02(kref,href);

 TimeIMBDF2=[runtime1,runtime2,runtime3,runtime4];
 [convIMBDF2,errorIMBDF2]=conv_Brusselator(soln1,soln2,soln3,soln4,solnref2,solnref3,solnref4,solnref5,step);
 
 %% results for IMEX-TR
[runtime1,soln1] =  Brusselator2D_IMEX_TR(k1,h1);
[runtime2,soln2] =  Brusselator2D_IMEX_TR(k2,h2);
[runtime3,soln3] = Brusselator2D_IMEX_TR(k3,h3);
[runtime4,soln4] = Brusselator2D_IMEX_TR(k4,h4);
%[~,solnref] = Brusselator2D_ETDpade02(kref,href);

 TimeIMTR=[runtime1,runtime2,runtime3,runtime4];
 [convIMTR,errorIMTR]=conv_Brusselator(soln1,soln2,soln3,soln4,solnref2,solnref3,solnref4,solnref5,step);
 
%% results for IMEX-Adams2
[runtime1,soln1] =  Brusselator2D_IMEX_Adams2(k1,h1);
[runtime2,soln2] =  Brusselator2D_IMEX_Adams2(k2,h2);
[runtime3,soln3] = Brusselator2D_IMEX_Adams2(k3,h3);
[runtime4,soln4] = Brusselator2D_IMEX_Adams2(k4,h4);
%[~,solnref] = Brusselator2D_ETDpade02(kref,href);

 TimeIMAD2=[runtime1,runtime2,runtime3,runtime4];
 [convIMAD2,errorIMAD2]=conv_Brusselator(soln1,soln2,soln3,soln4,solnref2,solnref3,solnref4,solnref5,step); 
 
 %% Effiency plot
 % Comparison for ETD Schemes
%  Time_mat = [TimeCN;TimeRRP;Timepa;TimeBDF2];
%  Error_mat = [errorCN;errorRRP;errorpa;errorBDF2];
 
 % Comparison for IMEX Schemes
 Time_mat = [TimeRDP;TimeIMBDF2;TimeIMTR;TimeIMAD2;TimeRDPS];
 Error_mat = [errorRDP;errorIMBDF2;errorIMTR;errorIMAD2;errorRDPS];
 
 
 efficiency_plot_Brusselator2D(Time_mat,Error_mat)
 convergence_plot_Brusselator2D(step,Error_mat)
 
 save BrusselatorDATA Time_mat Error_mat

%% Display results
fprintf(' Results for ETDRDP \n')
fprintf('k             h            error        conv       Time \n');
for i = 1:4
    fprintf('%.6f   %.6f    %1.4e      %.2f      %.5f\n', step(i),space(i),errorRDP(i),convRDP(i),TimeRDP(i))
end

fprintf(' Results for ETDRDP Split \n')
fprintf('k             h            error        conv       Time \n');
for i = 1:4
    fprintf('%.6f   %.6f    %1.4e      %.2f      %.5f\n', step(i),space(i),errorRDPS(i),convRDPS(i),TimeRDPS(i))
end
% fprintf('\n Results for ETDCN \n')
% fprintf('k             h            error        conv       Time \n');
% for i = 1:4
%     fprintf('%.6f   %.6f    %1.4e      %.2f      %.5f\n', step(i),space(i),errorCN(i),convCN(i),TimeCN(i))
% end
% 
% fprintf('\nResults for ETDpade02\n')
% fprintf('k             h            error        conv       Time\n');
% for i = 1:4
%     fprintf('%.6f   %.6f    %1.4e      %.2f      %.5f\n', step(i),space(i),errorpa(i),convpa(i),Timepa(i))
% end
% fprintf('\nResults for BDF2\n')
% fprintf('k             h            error        conv       Time\n');
% for i = 1:4
%     fprintf('%.6f   %.6f    %1.4e      %.2f      %.5f\n', step(i),space(i),errorBDF2(i),convBDF2(i),TimeBDF2(i))
% end
fprintf('\nResults for IMEX-BDF2\n')
fprintf('k             h            error        conv       Time\n');
for i = 1:4
    fprintf('%.6f   %.6f    %1.4e      %.2f      %.5f\n', step(i),space(i),errorIMBDF2(i),convIMBDF2(i),TimeIMBDF2(i))
end

fprintf('\nResults for IMEX-TR\n')
fprintf('k             h            error        conv       Time\n');
for i = 1:4
    fprintf('%.6f   %.6f    %1.4e      %.2f      %.5f\n', step(i),space(i),errorIMTR(i),convIMTR(i),TimeIMTR(i))
end

fprintf('\nResults for IMEX-Adams2\n')
fprintf('k             h            error        conv       Time\n');
for i = 1:4
    fprintf('%.6f   %.6f    %1.4e      %.2f      %.5f\n', step(i),space(i),errorIMAD2(i),convIMAD2(i),TimeIMAD2(i))
end




