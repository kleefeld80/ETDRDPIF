% function to compute error and convergence rates

function [conv,error]=conv_ekinetics(soln1,soln2,soln3,soln4,solnref2,solnref3,solnref4,solnref5,step)

 
%initialize variables
error = zeros(1,3); conv = error;
n1 = sqrt(length(soln1)); n2 = sqrt(length(soln2));n3 = sqrt(length(soln3)); 
n4 = sqrt(length(soln4)); n5 = sqrt(length(solnref5));m1 = n1^2;m2 = n2^2;
m3 = n3^2;m4 = n4^2;

%% Compute Error using fine grid solutions
% Error 1
% extract course grid solution from fine grid
% SOLTEMP = zeros(n1);
% SOL = reshape(solnref2,n2,n2);
% SOL = SOL'; 
% for i = 1:n1
%     for j = 1:n1
%        SOLTEMP(i,j) = SOL(2*i,2*j);
%     end
% end
%SOLTEMP = SOLTEMP';
%soln11 = reshape(SOLTEMP,m1,1);
soln11 = abs(soln1-solnref2);
error(1) = max(soln11);
sum1=0;
sum2 = 0;
for i = 1:m1
    sum1 = sum1 + soln11(i)^2;
    sum2 = sum2 + solnref2(i).^2;    
end
error(1) = sqrt(((1/(n1+1))^2)*sum1)/sqrt(((1/(n1+1))^2)*sum2);

% Error 2
% extract course grid solution from fine grid
% SOLTEMP = zeros(n2);
% SOL = reshape(solnref3,n3,n3);
% SOL = SOL'; 
% for i = 1:n2
%     for j = 1:n2
%        SOLTEMP(i,j) = SOL(2*i,2*j);
%     end
% end
% SOLTEMP = SOLTEMP';
% soln22 = reshape(SOLTEMP,m2,1);
soln22 = abs(soln2-solnref3);
%error(2) = max(soln22);

sum1=0;
sum2=0;
for i = 1:m2
    sum1 = sum1 + soln22(i)^2;
    sum2 = sum2 + solnref3(i)^2;
end
error(2) = sqrt(((1/(n2+1))^2)*sum1)/sqrt(((1/(n2+1))^2)*sum2);

% Error 3
% extract course grid solution from fine grid
% SOLTEMP = zeros(n3);
% SOL = reshape(solnref4,n4,n4);
% SOL = SOL'; 
% for i = 1:n3
%     for j = 1:n3
%        SOLTEMP(i,j) = SOL(2*i,2*j);
%     end
% end
% SOLTEMP = SOLTEMP';
% soln33 = reshape(SOLTEMP,m3,1);
soln33 = abs(soln3-solnref4);
%error(3) = max(soln33);

sum1=0;
sum2=0;
for i = 1:m3
    sum1 = sum1 + soln33(i)^2;
    sum2 = sum2 + solnref4(i)^2;
end
error(3) = sqrt(((1/(n3+1))^2)*sum1)/sqrt(((1/(n3+1))^2)*sum2);

%Error 4
% extract course grid solution from fine grid
% SOLTEMP = zeros(n4);
% SOL = reshape(solnref5,n5,n5);
% SOL = SOL'; 
% for i = 1:n4
%     for j = 1:n4
%        SOLTEMP(i,j) = SOL(2*i,2*j);
%     end
% end
% SOLTEMP = SOLTEMP';
% soln44 = reshape(SOLTEMP,m4,1);
soln44 = abs(soln4-solnref5);
%error(4) = max(soln44);

sum1=0;
sum2=0;
for i = 1:m4
    sum1 = sum1 + soln44(i)^2;
    sum2 = sum2 + solnref5(i)^2;
end
error(4) = sqrt(((1/(n4+1))^2)*sum1)/sqrt(((1/(n4+1))^2)*sum2);

%calculate rate of convergence
conv(1)=0;
for i = 2:4
    conv(i) = (log(error(i-1)/error(i)))/(log(step(i-1)/step(i)));
end
