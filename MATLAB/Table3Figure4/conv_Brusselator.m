% function to compute error and convergence rates

function [conv,error]=conv_Brusselator(soln1,soln2,soln3,soln4,solnref2,solnref3,solnref4,solnref5,step)

%initialize variables
error = zeros(1,4); conv = error;
n1 = sqrt(length(soln1)); n2 = sqrt(length(soln2));n3 = sqrt(length(soln3)); 
n4 = sqrt(length(soln4)); n5 = sqrt(length(solnref5));

%% Compute Error using fine grid solutions

% Error 1
%sol_temp = finetocoarse(n1,n2,2,solnref2);
sol_temp= solnref2;
sol_error = abs(sol_temp-soln1);
error(1) = max(sol_error);
% sum=0;
% for i = 1:n1
%     sum = sum + sol_error(i)^2;
% end
% error(1) = sqrt((1/(n1-1))*sum);

% Error 2
%sol_temp = finetocoarse(n2,n3,2,solnref3);
sol_temp= solnref3;
sol_error = abs(sol_temp-soln2);
error(2) = max(sol_error);
% sum=0;
% for i = 1:n2
%     sum = sum + sol_error(i)^2;
% end
% error(2) = sqrt((1/(n2-1))*sum);

% Error 3
%sol_temp = finetocoarse(n3,n4,2,solnref4);
sol_temp= solnref4;
sol_error = abs(sol_temp-soln3);
error(3) = max(sol_error);
% sum=0;
% for i = 1:n3
%     sum = sum + sol_error(i)^2;
% end
% error(3) = sqrt((1/(n3-1))*sum);

%Error 4
%sol_temp = finetocoarse(n4,n5,2,solnref5);
sol_temp= solnref5;
sol_error = abs(sol_temp-soln4);
error(4) = max(sol_error);
% sum=0;
% for i = 1:n4
%     sum = sum + sol_error(i)^2;
% end
% error(4) = sqrt((1/(n4-1))*sum);

% %Error 5
% sol_temp = zeros(n5,1);
% for i = 1:n5
%    sol_temp(i) = solnref(8*i);  
% end
% sol_error = abs(sol_temp-soln5);
% sum=0;
% for i = 1:n5
%     sum = sum + sol_error(i)^2;
% end
% error(5) = sqrt(((1/(n5+1)))*sum);

%calculate rate of convergence
conv(1)=0;
for i = 2:4
    conv(i) = (log(error(i-1)/error(i)))/(log(step(i-1)/step(i)));
end
