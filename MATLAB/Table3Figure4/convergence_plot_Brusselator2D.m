% Code to produce convergence plot 

function convergence_plot_Brusselator2D(step,Error_mat)
% INPUT: Time_mat is matrix whose rows correspond to the time values of
% each algorithm
%        Error_mat is a matrix whose rows correspond to the error values of
%        each algorithm. Must be the same size as Time_mat

%% Extract Data

% Extract Error
 error1= Error_mat(1,:);
 error2= Error_mat(2,:);
 error3= Error_mat(3,:);
 error4= Error_mat(4,:);
 error5= Error_mat(5,:);
%  error6= Error_mat(6,:);
% error7= Error_mat(7,:);


%% Produce Plots
mk = {'o','s','^','h','v','+','d','<','*','>'};
c = [0 0 0;0 0 1;1 0 0;1 0 1;.75 0.25 0.25;0.25 0.5 0;0.75 0.5 0;0.75 0.25 1;1 0.25 0.75];
figure
set(0,'defaultaxesfontsize',14)
loglog(step,error1,'s-','color',c(1,:),...
    'LineWidth',2.5,...
    'MarkerSize',10')
hold on
loglog(step,error2,'o--','color',c(2,:),...
    'LineWidth',2.5,...
    'MarkerSize',10')
loglog(step,error3,'^--','color',c(3,:),...
    'LineWidth',2.5,...
    'MarkerSize',10')
loglog(step,error4,'h--','color',c(4,:),...
    'LineWidth',2.5,...
    'MarkerSize',10')
loglog(step,error5,'v--','color',c(6,:),...
    'LineWidth',2.5,...
    'MarkerSize',10')
% loglog(step,error6,'<--','color',c(7,:),...
%     'LineWidth',2.5,...
%     'MarkerSize',10')
% semilogy(step,error7,'--ms',...
%     'LineWidth',3,...
%     'MarkerSize',5,...
%     'MarkerEdgeColor','k',...
%     'MarkerFaceColor','m')
grid on
%legend('ETD-RDP','ETD-CN','ETD-P02')
%legend('ETD-RDP','ETD-CN','ETD-P02','IMEX-BDF2')
%legend('ETDRDP','BDF2','ROB2','IRKLA','SDIRK','IMEX-BDF2')
legend('ETD-RDP','IMEX-BDF2','IMEX-TR','IMEX-Adam2','ETD-RDPsplit')



xlabel('\bf\fontsize{14} Time Step')
ylabel('\bf\fontsize{14} Error')




