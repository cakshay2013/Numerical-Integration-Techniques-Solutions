% Let A0 = [A011 A012 A013 A014   and K = [k^2     0  0   sqrt(k)
%           A021 A022 A023 A024            0       0  0   0
%           A031 A032 A033 A034            0       0  0   0
%           A041 A042 A043 A044]           sqrt(k) 0  0   k]
% where k=2

% As A=A0*(eye(4)+K), A= A0+A0*K

% A = [A011+A011*k^2+A014*sqrt(k)  A012   A013   A014+A011*sqrt(k)+A014*k
%      A021+A021*k^2+A024*sqrt(k)  A022   A023   A024+A021*sqrt(k)+A024*k
%      A031+A031*k^2+A034*sqrt(k)  A032   A033   A034+A031*sqrt(k)+A034*k
%      A041+A041*k^2+A044*sqrt(k)  A042   A043   A044+A041*sqrt(k)+A044*k]

% A' = [A011*k*2+A014*(1/(2*sqrt(k)))  0   0  A011*(1/(2*sqrt(k)))+A014
%       A021*k*2+A024*(1/(2*sqrt(k)))  0   0  A021*(1/(2*sqrt(k)))+A024
%       A031*k*2+A034*(1/(2*sqrt(k)))  0   0  A031*(1/(2*sqrt(k)))+A034
%       A041*k*2+A044*(1/(2*sqrt(k)))  0   0  A041*(1/(2*sqrt(k)))+A044]
% where A'= dA/dk

% let K1= [2*k            0   0   1/(2*sqrt(k))
%          0              0   0   0
%          0              0   0   0
%          1/(2*sqrt(k))  0   0   1]
% where K1 = dK/dk

% A' or dA/dk can be written as, dA/dk=P1=A0*K1

% Eigenvalues of A matrix are all real. Thus modes are anologous to
% eigenvalues. Thus, the eigenvalue sensitivities gives us sensitivities of
% modes w.r.t. the variable k.
% columns of participation matrix P of A gives us idea of modes

% Thus modal sensitivities or dmode(a)/dk = ΣΣPHI(j,a)*PSI(a,i)*P1(i,j)
% where i,i varies from 1 to n and n=4.
% storing dmode(a)/dk values in an array Z(a) where a=1 to n.

% by sorting in descending order the absolute values of Z, we get the mode that is most
% sensitive.

% thus, [B,I] = sort(abs(Z),'descend'), where I is 1*n array of sorted
% indices of Z

% P(:,I(1)) gives us the required answer, where P is the participation
% matrix

%***************************************************************************
clc;
data=xlsread("EE5230_2021_Quiz3_DS.xlsx", 17,'B9:B24');
A0=zeros(4,4);
ind=1;
for k=1:4
    for l=1:4
        A0(k,l)=data(ind);
        ind=ind+1;
    end
end
k=2;
K = [k^2 0 0 sqrt(k); 0 0 0 0; 0 0 0 0; sqrt(k) 0 0 k];
A = A0*(eye(4)+K);
% dA/dk= A0*dK/dk
n=4;
[PHI,Diag] = eig(A);
eig_val = zeros(n,1);
for i = 1:n
    eig_val(i) = Diag(i,i);
end
eig_val
PHI
PSI = inv(PHI)
P = PHI.*PSI % Participation factor matrix
K1 = [2*k 0 0 1/(2*sqrt(k)); 0 0 0 0; 0 0 0 0; 1/(2*sqrt(k)) 0 0 1]; % dK/dk
P1=A0*K1; % dA/dk= A0*dK/dk
Z=zeros(1,4);
% dlambda(k)/dk or eigenvalue/mode sensitivities:
% dlambda(k)/dk = dA(i,j)/dk*PHI(j,k)*PSI(k,i) or P1(i,j)*PHI(j,a)*PSI(a,i)
for i=1:n
    for j=1:n
        for a=1:n
            Z(a)= P1(i,j)*PHI(j,a)*PSI(a,i);
        end
    end
end
Z;
[B,I] = sort(abs(Z),'descend');
B;
I;
P(:,I(1))
% row of P gives states info and column give mode/eigenvalue info