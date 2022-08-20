% Eulers, Modified Eulers and RK4 Solution for linear/non-linear differential-algebriac equations of the form:
% dx/dt=f(x,z), g(x,z)=0
clc;
n = 3; % n is the size of x,z vectors for each time instant. n = iterations+1 = final_time/h + 1
h = 0.2; % h is the step size

x1(1,1) = 2.52;
x2(1,1) = 4.09;
z1(1,1) = 1.42;

b1 = x1(1,1);
b2 = x2(1,1);
b3 = z1(1,1);

c = (z1(1,1)^5)*sin(z1(1,1))+1.5*(z1(1,1)^3)+x1(1,1)*x2(1,1)*z1(1,1)+x1(1,1)^2;
%x1(1,1)*(z1(1,1)^5) + x2(1,1)*(z1(1,1)^3) + x1(1,1)*x2(1,1)*z1(1,1);

[a1,a2,a3] = Euler_DAE(b1,b2,b3,h,n,c);
[a4,a5,a6] = MEuler_DAE(b1,b2,b3,h,n,c);
[a7,a8,a9] = RK4_DAE(b1,b2,b3,h,n,c);

[a1,a2,a3].'
[a4,a5,a6].'
[a7,a8,a9].'

function [a1,a2,a3] = Euler_DAE(b1,b2,b3,h,n,c)
x1 = zeros(n,1);
x2 = zeros(n,1);
%x3 = zeros(n,1);
z1 = zeros(n,1);
%z2 = zeros(n,1);
x1(1,1) = b1;
x2(1,1) = b2;
z1(1,1) = b3;
for i = 1:(n-1)
    x1(i+1,1) = x1(i,1) + h*F1(x1(i,1),x2(i,1),NR_DAE(z1(i,1),x1(i,1),x2(i,1),c));
    x2(i+1,1) = x2(i,1) + h*F2(x1(i,1),x2(i,1),NR_DAE(z1(i,1),x1(i,1),x2(i,1),c));
    z1(i+1,1) = NR_DAE(z1(i,1),x1(i+1,1),x2(i+1,1),c);
end
a1 = x1;
a2 = x2;
a3 = z1;
end

function [a4,a5,a6] = MEuler_DAE(b1,b2,b3,h,n,c)
m = 2; % m is the iterations in the inner corrector loop
x1 = zeros(n,1);
x2 = zeros(n,1);
%x3 = zeros(n,1);
z1 = zeros(n,1);
%z2 = zeros(n,1);
xp1 = zeros(m,1);
xp2 = zeros(m,1);
zp1 = zeros(m,1);
x1(1,1) = b1;
x2(1,1) = b2;
z1(1,1) = b3;
for i = 1:(n-1)
    xp1(1,1) = x1(i,1) + h*F1(x1(i,1),x2(i,1),NR_DAE(z1(i,1),x1(i,1),x2(i,1),c));
    xp2(1,1) = x2(i,1) + h*F2(x1(i,1),x2(i,1),NR_DAE(z1(i,1),x1(i,1),x2(i,1),c));
    zp1(1,1) = NR_DAE(z1(i,1),xp1(1,1),xp2(1,1),c);
for j = 1:(m-1)
if m == 1
    break
end
    xp1(j+1,1) = x1(i,1) + 0.5*h*(F1(x1(i,1),x2(i,1),NR_DAE(z1(i,1),x1(i,1),x2(i,1),c)) + F1(xp1(j,1),xp2(j,1),NR_DAE(zp1(j,1),xp1(j,1),xp2(j,1),c)));
    xp2(j+1,1) = x2(i,1) + 0.5*h*(F2(x1(i,1),x2(i,1),NR_DAE(z1(i,1),x1(i,1),x2(i,1),c)) + F2(xp1(j,1),xp2(j,1),NR_DAE(zp1(j,1),xp1(j,1),xp2(j,1),c)));
    zp1(j+1,1) = NR_DAE(zp1(j,1),xp1(j+1,1),xp2(j+1,1),c);
end
    x1(i+1,1) = x1(i,1) + 0.5*h*(F1(x1(i,1),x2(i,1),z1(i,1)) + F1(xp1(m,1),xp2(m,1),NR_DAE(zp1(m,1),xp1(m,1),xp2(m,1),c)));
    x2(i+1,1) = x2(i,1) + 0.5*h*(F2(x1(i,1),x2(i,1),z1(i,1)) + F2(xp1(m,1),xp2(m,1),NR_DAE(zp1(m,1),xp1(m,1),xp2(m,1),c)));
    z1(i+1,1) = NR_DAE(zp1(m,1),x1(i+1,1),x2(i+1,1),c);
end
a4 = x1;
a5 = x2;
a6 = z1;
end

function [a7,a8,a9] = RK4_DAE(b1,b2,b3,h,n,c)
x1 = zeros(n,1);
x2 = zeros(n,1);
%x3 = zeros(n,1);
z1 = zeros(n,1);
%z2 = zeros(n,1);
c1 = zeros(2,n-1);
c2 = zeros(2,n-1);
c3 = zeros(2,n-1);
c4 = zeros(2,n-1);
x1(1,1) = b1;
x2(1,1) = b2;
z1(1,1) = b3;
for i = 1:(n-1)
    c1(1,i) = h*F1(x1(i,1),x2(i,1),NR_DAE(z1(i,1),x1(i,1),x2(i,1),c));
    c1(2,i) = h*F2(x1(i,1),x2(i,1),NR_DAE(z1(i,1),x1(i,1),x2(i,1),c));
    
    c2(1,i) = h*F1(x1(i,1)+0.5*c1(1,i),x2(i,1)+0.5*c1(2,i),NR_DAE(z1(i,1),x1(i,1)+0.5*c1(1,i),x2(i,1)+0.5*c1(2,i),c));
    c2(2,i) = h*F2(x1(i,1)+0.5*c1(1,i),x2(i,1)+0.5*c1(2,i),NR_DAE(z1(i,1),x1(i,1)+0.5*c1(1,i),x2(i,1)+0.5*c1(2,i),c));
    
    c3(1,i) = h*F1(x1(i,1)+0.5*c2(1,i),x2(i,1)+0.5*c2(2,i),NR_DAE(z1(i,1),x1(i,1)+0.5*c2(1,i),x2(i,1)+0.5*c2(2,i),c));
    c3(2,i) = h*F2(x1(i,1)+0.5*c2(1,i),x2(i,1)+0.5*c2(2,i),NR_DAE(z1(i,1),x1(i,1)+0.5*c2(1,i),x2(i,1)+0.5*c2(2,i),c));

    c4(1,i) = h*F1(x1(i,1)+c3(1,i),x2(i,1)+c3(2,i),NR_DAE(z1(i,1),x1(i,1)+c3(1,i),x2(i,1)+c3(2,i),c));
    c4(2,i) = h*F2(x1(i,1)+c3(1,i),x2(i,1)+c3(2,i),NR_DAE(z1(i,1),x1(i,1)+c3(1,i),x2(i,1)+c3(2,i),c));

    x1(i+1,1) = x1(i,1) + (1/6)*(c1(1,i)+2*c2(1,i)+2*c3(1,i)+c4(1,i));
    x2(i+1,1) = x2(i,1) + (1/6)*(c1(2,i)+2*c2(2,i)+2*c3(2,i)+c4(2,i));
    z1(i+1,1) = NR_DAE(z1(i,1),x1(i+1,1),x2(i+1,1),c);
end
a7 = x1;
a8 = x2;
a9 = z1;
end

function z = NR_DAE(Z_1,x1,x2,c)
for j = 1:2 % 2 NR iterations
    Z_1 = Z_1 - G1(x1,x2,Z_1,c)/DG1(x1,x2,Z_1);
end
z = Z_1;
end

function f1 = F1(x1,x2,z1) %f(x(i), z(i))
f1 = -x1*(z1^2)-x1*x2;
%-x1^2-x2+z1;
end

function f2 = F2(x1,x2,z1) %f(x(i), z(i))
f2 = -x1+(x2^2)*z1;
%-x1+3*x2*z1; 
end

function g1 = G1(x1,x2,z1,c) %g(x(i), z(i))
g1 = (z1^5)*sin(z1)+1.5*(z1^3)+x1*x2*z1+x1^2-c;
%x1*(z1^5)+x2*(z1^3)+x1*x2*z1-c;
end

function dg1 = DG1(x1,x2,z1) %dg(x(i), z(i))/dz
dg1 = 5*(z1^4)*sin(z1)+(z1^5)*cos(z1)+3*1.5*(z1^2)+x1*x2;
%5*x1*(z1^4)+3*x2*(z1^2)+x1*x2;
end