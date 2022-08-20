% Eulers, Modified Eulers and RK4 Solution for linear/non-linear differential equations of the form:
% dx/dt=f(x,z)
clc;
n = 3; % n is the size of x,z vectors for each time instant. n = iterations+1 = final_time/h + 1
h = 0.5; % h is the step size

x1(1,1) = 1;
x2(1,1) = 0.5;

b1 = x1(1,1);
b2 = x2(1,1);
%x1(1,1)*(z1(1)^5) + x2(1,1)*(z1(1)^3) + x1(1,1)*x2(1,1)*z1(1);

[a1,a2] = Euler_DAE(b1,b2,h,n);
[a3,a4] = MEuler_DAE(b1,b2,h,n);
[a5,a6] = RK4_DAE(b1,b2,h,n);

[a1,a2].'
[a3,a4].'
[a5,a6].'

function [a1,a2] = Euler_DAE(b1,b2,h,n)
x1 = zeros(n,1);
x2 = zeros(n,1);
%x3 = zeros(n,1);
x1(1,1) = b1;
x2(1,1) = b2;
for i = 1:(n-1)
    x1(i+1,1) = x1(i,1) + h*F1(x1(i,1),x2(i,1));
    x2(i+1,1) = x2(i,1) + h*F2(x1(i,1),x2(i,1));
end
a1 = x1;
a2 = x2;
end

function [a3,a4] = MEuler_DAE(b1,b2,h,n)
x1 = zeros(n,1);
x2 = zeros(n,1);
%x3 = zeros(n,1);
x1(1,1) = b1;
x2(1,1) = b2;
m = 2; % m is the iterations in the inner corrector loop
xp1 = zeros(m,1);
xp2 = zeros(m,1);
for i = 1:(n-1)
    xp1(1,1) = x1(i,1) + h*F1(x1(i,1),x2(i,1));
    xp2(1,1) = x2(i,1) + h*F2(x1(i,1),x2(i,1));
for j = 1:(m-1)
if m == 1
    break
end
    xp1(j+1,1) = x1(i,1) + 0.5*h*(F1(x1(i,1),x2(i,1)) + F1(xp1(j,1),xp2(j,1)));
    xp2(j+1,1) = x2(i,1) + 0.5*h*(F2(x1(i,1),x2(i,1)) + F2(xp1(j,1),xp2(j,1)));
end
    x1(i+1,1) = x1(i,1) + 0.5*h*(F1(x1(i,1),x2(i,1)) + F1(xp1(m,1),xp2(m,1)));
    x2(i+1,1) = x2(i,1) + 0.5*h*(F2(x1(i,1),x2(i,1)) + F2(xp1(m,1),xp2(m,1)));
end
a3 = x1;
a4 = x2;
end

function [a5,a6] = RK4_DAE(b1,b2,h,n)
x1 = zeros(n,1);
x2 = zeros(n,1);
%x3 = zeros(n,1);
x1(1,1) = b1;
x2(1,1) = b2;

c1 = zeros(2,n-1);
c2 = zeros(2,n-1);
c3 = zeros(2,n-1);
c4 = zeros(2,n-1);

for i = 1:(n-1)
    c1(1,i) = h*F1(x1(i,1),x2(i,1));
    c1(2,i) = h*F2(x1(i,1),x2(i,1));
    
    c2(1,i) = h*F1(x1(i,1)+0.5*c1(1,i),x2(i,1)+0.5*c1(2,i));
    c2(2,i) = h*F2(x1(i,1)+0.5*c1(1,i),x2(i,1)+0.5*c1(2,i));
    
    c3(1,i) = h*F1(x1(i,1)+0.5*c2(1,i),x2(i,1)+0.5*c2(2,i));
    c3(2,i) = h*F2(x1(i,1)+0.5*c2(1,i),x2(i,1)+0.5*c2(2,i));

    c4(1,i) = h*F1(x1(i,1)+c3(1,i),x2(i,1)+c3(2,i));
    c4(2,i) = h*F2(x1(i,1)+c3(1,i),x2(i,1)+c3(2,i));

    x1(i+1,1) = x1(i,1) + (1/6)*(c1(1,i)+2*c2(1,i)+2*c3(1,i)+c4(1,i));
    x2(i+1,1) = x2(i,1) + (1/6)*(c1(2,i)+2*c2(2,i)+2*c3(2,i)+c4(2,i));
end
a5 = x1;
a6 = x2;
end

function f1 = F1(x1,x2) %f(x(i), z(i))
f1 = x1-x2*x1;
%-x1^2-x2+z1;
end

function f2 = F2(x1,x2) %f(x(i), z(i))
f2 = x1*x1 - x2;
%-x1+3*x2*z1; 
end