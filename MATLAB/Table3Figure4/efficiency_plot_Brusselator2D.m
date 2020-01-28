% Code to produce efficiency plot algorithms

function efficiency_plot_Brusselator2D(Time_mat,Error_mat)
% INPUT: Time_mat is matrix whose rows correspond to the time values of
% each algorithm
%        Error_mat is a matrix whose rows correspond to the error values of
%        each algorithm. Must be the same size as Time_mat

%% Extract Data
% Extract Time
time1= Time_mat(1,:);
time2= Time_mat(2,:);
time3= Time_mat(3,:);
time4= Time_mat(4,:);
time5= Time_mat(5,:);
%time6= Time_mat(6,:);



% Extract Error
error1= Error_mat(1,:);
error2= Error_mat(2,:);
error3= Error_mat(3,:);
error4= Error_mat(4,:);
error5= Error_mat(5,:);
% error6= Error_mat(6,:);

siz=8;
ls=2;
%% Produce Plots
mk = {'o','s','^','h','v','+','d','<','*','>'};
c = [0 0 0;0 0 1;1 0 0;.25 .5 .25;.75 0.25 0.25;0.25 0.5 0;0.75 0.5 0;0.75 0.25 1;1 0.25 0.75];
set(0,'defaultaxesfontsize',12)
loglog(error1,time1,'s-','color',c(1,:),...
    'LineWidth',ls,...
    'MarkerSize',siz)
hold on
loglog(error2,time2,'o-.','color',c(2,:),...
    'LineWidth',ls,...
    'MarkerSize',siz)
loglog(error3,time3,'p-.','color',c(3,:),...
    'LineWidth',ls,...
    'MarkerSize',siz)
loglog(error4,time4,'<--','color',c(5,:),...
    'LineWidth',ls,...
    'MarkerSize',siz)
loglog(error5,time5,'v-.','color',c(7,:),...
    'LineWidth',ls,...
    'MarkerSize',siz)
% loglog(error6,time6,'h-','color',c(7,:),...
%     'LineWidth',2.5,...
%     'MarkerSize',10)
grid on
%legend('ETDRDP','BDF2','ROB2','IRKLA','SDIRK','IMEX-BDF2')
%legend('ETD-RDP','ETD-CN','ETD-P02')
%legend('ETD-RDP','ETD-CN','ETD-P02','IMEX-BDF2')
legend('ETD-RDP','IMEX-BDF2','IMEX-TR','IMEX-Adam2','ETD-RDP-IF')

xlabel('error')
ylabel('CPU time in seconds')
%xlabel('\bf\fontsize{14} Time Step')
%ylabel('\bf\fontsize{14} Error')
set(gca,'LineWidth', 1);
set(gca,'FontSize',10);
set(gca,'FontWeight','bold');
axis([1e-3,1e-0,1e-2,1e1])
print -depsc2 efficiency2.eps

%xlabel('\bf\fontsize{12} Error')
%ylabel('\bf\fontsize{12}CPU Time')




