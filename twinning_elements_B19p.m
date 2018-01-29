function [n,a,eta,K,shear] = twinning_elements_B19p(Uj,Ui,kappa,f1m,f2m,f3m)
% C
C = inv(Uj)*(Ui*Ui)*inv(Uj);

% eigenvalues and eigenvectors of C in ascending order
[evec,eval]=eig(C);
temp=[transpose(evec(:,1)) eval(1,1);transpose(evec(:,2)) eval(2,2);transpose(evec(:,3)) eval(3,3)];
[temp1,temp2]=sort(temp(:,4));
temp_new=temp(temp2,:);
lam1=temp_new(1,4);  lam2=temp_new(2,4);  lam3=temp_new(3,4);
e1=transpose(temp_new(1,1:3));  e2=transpose(temp_new(2,1:3));  e3=transpose(temp_new(3,1:3));

% rho
rho=norm((sqrt(lam3)-sqrt(lam1))/sqrt(lam3-lam1)*(-sqrt(1-lam1)*Uj*e1+kappa*sqrt(lam3-1)*Uj*e3));


%% a
a=rho*(sqrt(lam3*(1-lam1)/(lam3-lam1))*e1+kappa*sqrt(lam1*(lam3-1)/(lam3-lam1))*e3);


%% eta
eta_temp=inv(Uj)*unit(a);
temp=max(abs(eta_temp));
eta=eta_temp/temp;


%% n
n=unit((sqrt(lam3)-sqrt(lam1))/sqrt(lam3-lam1)*(-sqrt(1-lam1)*Uj*e1+kappa*sqrt(lam3-1)*Uj*e3));  % B2



%% K
l1=[1;0;0];  l2=[0;1;0];  l3=[0;0;1];
[l1_r,l2_r,l3_r]=reciprocal(l1,l2,l3);
K_temp=unit(n);

% 
Rec = [l1_r l2_r l3_r];
hkl = eye(3)/Rec*K_temp;
K=[hkl(1) hkl(2) hkl(3)];
temp=max(abs(K));
K=K/temp;

I=eye(3);
[f1m_r,f2m_r,f3m_r]=reciprocal(f1m,f2m,f3m);
K_temp=unit(I/Uj*n);
% syms h k l
% [h k l]=solve(h*f1m_r(1)+k*f2m_r(1)+l*f3m_r(1)==K_temp(1),h*f1m_r(2)+k*f2m_r(2)+l*f3m_r(2)==K_temp(2),h*f1m_r(3)+k*f2m_r(3)+l*f3m_r(3)==K_temp(3),h,k,l);
% K=double([h;k;l]);

hkl = (eye(3) / [f1m_r f2m_r f3m_r]) * K_temp;
K=[hkl(1);hkl(2);hkl(3)];
% temp=max(abs(K));
% K=K/temp; 
n=K;  % ***********



% syms h k l
% [h k l]=solve(h*f1m(1)+k*f2m(1)+l*f3m(1)==a(1),h*f1m(2)+k*f2m(2)+l*f3m(2)==a(2),h*f1m(3)+k*f2m(3)+l*f3m(3)==a(3),h,k,l);
% ee=double([h;k;l]);
hkl = (eye(3) / [f1m f2m f3m]) * a;
ee = [hkl(1);hkl(2);hkl(3)];
% temp=max(abs(K));
% K=K/temp; 
a=ee;  % ***********


%% shear
shear=sqrt(lam3)-sqrt(lam1);