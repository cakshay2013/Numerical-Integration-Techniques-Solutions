% Solution for time scale decomposition of linear differential equations of the form:
% dx/dt=f(x)

x01=4.5;
x02=4.5;
x03=4.5;
x04=4.5;
A = [-3.39 -5.48 -1.31 -5.3;-5.47 -2.59 -6.09 -1.83;-18.1 -25.26 -18.5 -34.61;-38.92 -16.81 -29.13 -23.08];
[C,d] = Decompose(x01,x02,x03,x04,A)

function [C,d] = Decompose(x01,x02,x03,x04,A)
p = zeros(4);
q = zeros(4);
B = zeros(4);

p(1) = abs(F1(x01,x02,x03,x04,A));
p(2) = abs(F2(x01,x02,x03,x04,A));
p(3) = abs(F3(x01,x02,x03,x04,A));
p(4) = abs(F4(x01,x02,x03,x04,A));

[B,I] = sort(p,'ascend'); % I is index vector, do not edit this
alpha = B(3)/B(4);

f=2; %no. of state variables to be taken as fast variables
% here 3,4 are the indices of fast variables and 1,2 for slow variables
A1 = A(I(1):I(2), I(1):(2)); % slow-slow submatrix
A2 = A(I(1):I(2), I(3):(4)); % slow-fast submatrix
A3 = A(I(3):I(4), I(1):(2)); % fast-slow submatrix
A4 = A(I(3):I(4), I(3):(4)); % fast-fast submatrix

q(1) = x01;
q(2) = x02;
q(3) = x03;
q(4) = x04;

xfast = [q(I(3)); q(I(4))];
H0 = -inv(A4)*A3;
H1 = -inv(A4)*inv(A4)*A3*(A1 + A2*H0);
H = H0 + H1;
C = A4;
d = A3*xfast + A4*H*xfast;
end

function f1 = F1(x1,x2,x3,x4,A) %f(x(i), z(i))
f1 = A(1,:)*[x1; x2; x3; x4];
%-3.39*x1-5.48*x2-1.31*x3-5.3*x4;
%-x1^2-x2+z1;
end

function f2 = F2(x1,x2,x3,x4,A) %f(x(i), z(i))
f2 = A(2,:)*[x1; x2; x3; x4];
%-5.47*x1-2.59*x2-6.09*x3-1.83*x4;
%-x1+3*x2*z1; 
end

function f3 = F3(x1,x2,x3,x4,A) %f(x(i), z(i))
f3 = A(3,:)*[x1; x2; x3; x4];
%-18.1*x1-25.26*x2-18.5*x3-34.61*x4;
%-x1+3*x2*z1; 
end

function f4 = F4(x1,x2,x3,x4,A) %f(x(i), z(i))
f4 = A(4,:)*[x1; x2; x3; x4];
%-38.92*x1-16.81*x2-29.13*x3-23.08*x4;
%-x1+3*x2*z1; 
end

%function f = F(x1,x2,x3,x4)
%A = [-3.39 -5.48 -1.31 -5.3;-5.47 -2.59 -6.09 -1.83;-18.1 -25.26 -18.5 -34.61;-38.92 -16.81 -29.13 -23.08];
%f = A*[x1; x2; x3; x4];
%end