% Solution for algebriac equations of the form:
% f(x,z)=0, g(x,z)=0

n = 3; % n is the size of x,z vectors for each time instant. n = iterations+1 = final_time/h + 1
x1 = zeros(n);
x2 = zeros(n);
%x3 = zeros(n);
z1 = zeros(n);
%z2 = zeros(n);

x1(1) = 1.8;
x2(1) = 4.2;
z1(1) = 0.9;
c = 0;
%x1(1)*(z1(1)^5) + x2(1)*(z1(1)^3) + x1(1)*x2(1)*z1(1);

[a1,a2,a3] = DAE_NL(x1,x2,z1,n);
[a1,a2,a3].'

function [a1,a2,a3] = DAE_NL(x1,x2,z1,n)
for i = 1:(n-1)
    J = j(x1(i),x2(i),z1(i));
    f = zeros(3,1);
    [f(1,1), f(2,1), f(3,1)] = FLOW(x1(i),x2(i),z1(i));
    K = inv(J)*f;
    x1(i+1) = x1(i) - K(1,:);
    x2(i+1) = x2(i) - K(2,:);
    z1(i+1) = z1(i) - K(3,:);
end
a1 = x1(n);
a2 = x2(n);
a3 = z1(n);
end

function [flow1, flow2, flow3] = FLOW(x1,x2,z1)
[flow1, flow2] = F(x1,x2,z1);
flow3 = G(x1,x2,z1);
end

function J = j(x1,x2,z1)
J = zeros(3,3);
[J(1,1), J(1,2), J(2,1), J(2,2)] = DFx(x1,x2,z1);
[J(1,3), J(2,3)] = DFz(x1,x2,z1);
[J(3,1), J(3,2)] = DGx(x1,x2,z1);
J(3,3) = DGz(x1,x2,z1);
end

function [f1, f2] = F(x1,x2,z1) %f(x(i), z(i))
f1 = x1+sqrt(x2)+x1*z1-6;
f2 = x1*x1+x2*z1-8;
%-x1^2-x2+z1;
end

function g = G(x1,x2,z1) %g(x(i), z(i))
g = x1 + x2*x1 - x1*z1 - 8;
%x1*(z1^5)+x2*(z1^3)+x1*x2*z1-c;
end

function [dfx1,dfx2,dfx3,dfx4] = DFx(x1,x2,z1) %dg(x(i), z(i))/dz
dfx1 = 1+z1;
dfx2 = 0.5/sqrt(x2);
dfx3 = 2*x1;
dfx4 = z1;
%2*x2*z1 + 0.5/sqrt(z1);
%5*x1*(z1^4)+3*x2*(z1^2)+x1*x2;
end

function [dfz1, dfz2] = DFz(x1,x2,z1) %dg(x(i), z(i))/dz
dfz1 = x1;
dfz2 = x2;
%2*x2*z1 + 0.5/sqrt(z1);
%5*x1*(z1^4)+3*x2*(z1^2)+x1*x2;
end

function [dgx1, dgx2] = DGx(x1,x2,z1) %dg(x(i), z(i))/dz
dgx1 = 1+x2-z1;
dgx2 = x1;
%2*x2*z1 + 0.5/sqrt(z1);
%5*x1*(z1^4)+3*x2*(z1^2)+x1*x2;
end

function dgz = DGz(x1,x2,z1) %dg(x(i), z(i))/dz
dgz = -x1;
%2*x2*z1 + 0.5/sqrt(z1);
%5*x1*(z1^4)+3*x2*(z1^2)+x1*x2;
end