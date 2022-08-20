x1 = 1.8;
x2 = 4.2;
z1 = 0.9;
c = 0;

j(x1,x2,z1)

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

function [dfx11,dfx12,dfx21,dfx22] = DFx(x1,x2,z1) %dg(x(i), z(i))/dz
dfx11 = 1+z1;
dfx12 = 0.5/sqrt(x2);
dfx21 = 2*x1;
dfx22 = z1;
%2*x2*z1 + 0.5/sqrt(z1);
%5*x1*(z1^4)+3*x2*(z1^2)+x1*x2;
end

function [dfz13, dfz23] = DFz(x1,x2,z1) %dg(x(i), z(i))/dz
dfz13 = x1;
dfz23 = x2;
%2*x2*z1 + 0.5/sqrt(z1);
%5*x1*(z1^4)+3*x2*(z1^2)+x1*x2;
end

function [dgx31, dgx32] = DGx(x1,x2,z1) %dg(x(i), z(i))/dz
dgx31 = 1+x2-z1;
dgx32 = x1;
%2*x2*z1 + 0.5/sqrt(z1);
%5*x1*(z1^4)+3*x2*(z1^2)+x1*x2;
end

function dgz33 = DGz(x1,x2,z1) %dg(x(i), z(i))/dz
dgz33 = -x1;
%2*x2*z1 + 0.5/sqrt(z1);
%5*x1*(z1^4)+3*x2*(z1^2)+x1*x2;
end