function [bplus,mplus,bminus,mminus,lambda,shapestrain]=habit_elements_HS(Uj,n,a,type)
%% m, b, and lamda
delta = dot(a,Uj*inv(Uj^2-eye(3))*n);  %display(delta)
lambda = 0.5*(1-sqrt(1+2/delta));
eta_new = trace(Uj^2)-det(Uj^2)-2+norm(a)^2/(2*delta);  %display(eta_new)


if delta<=-2 && eta_new>=0
    
    C=transpose(Uj+lambda*dyad(a,n))*(Uj+lambda*dyad(a,n));
    
    % eigenvalues and eigenvectors of C in ascending order
    [evec,eval]=eig(C);
    temp=[transpose(evec(:,1)) eval(1,1);transpose(evec(:,2)) eval(2,2);transpose(evec(:,3)) eval(3,3)];
    [temp1,temp2]=sort(temp(:,4));
    temp_new=temp(temp2,:);
    lam1=temp_new(1,4);  lam2=temp_new(2,4);  lam3=temp_new(3,4);
    e1=transpose(temp_new(1,1:3));  e2=transpose(temp_new(2,1:3));  e3=transpose(temp_new(3,1:3));
    
    
    %% mplus
    kappa=1;
    m_temp=(sqrt(lam3)-sqrt(lam1))/sqrt(lam3-lam1)*(-sqrt(1-lam1)*e1+kappa*sqrt(lam3-1)*e3);
    rho=norm(m_temp);
    mplus=m_temp/rho;
    
    
    %% bplus
    bplus=rho*(sqrt(lam3*(1-lam1)/(lam3-lam1))*e1+kappa*sqrt(lam1*(lam3-1)/(lam3-lam1))*e3);
    
    
    %% mminus
    kappa=-1;
    m_temp=(sqrt(lam3)-sqrt(lam1))/sqrt(lam3-lam1)*(-sqrt(1-lam1)*e1+kappa*sqrt(lam3-1)*e3);
    rho=norm(m_temp);
    mminus=m_temp/rho;
    
    
    %% bminus
    bminus=rho*(sqrt(lam3*(1-lam1)/(lam3-lam1))*e1+kappa*sqrt(lam1*(lam3-1)/(lam3-lam1))*e3);
    
    
    %% shapestrain
    shapestrain=norm(bplus);
    
    
else  % if A/M not possible
    bplus=zeros(3,1); mplus=bplus; bminus=bplus; mminus=bplus;
    lambda=0; shapestrain=0;
end