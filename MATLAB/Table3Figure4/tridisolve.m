% An implementation of Thomas algorithm
function x = tridisolve(A,d)


x = d;
n = length(x);
b = zeros(1,n);
a = zeros(1,n-1);
c = zeros(1,n-1);

for i = 1:n
    if i==n
     b(i) = A(i,i);
    else
    b(i) = A(i,i);
    a(i) = A(i+1,i);
    c(i) = A(i,i+1);
    end
end
    

%forward sweep
for j = 1:n-1
mu = a(j)/b(j);
b(j+1) = b(j+1) - mu*c(j);
x(j+1) = x(j+1) - mu*x(j);
end

%Backward Sweep
x(n) = x(n)/b(n);

for j = n-1:-1:1
x(j) = (x(j)-c(j)*x(j+1))/b(j);
end
