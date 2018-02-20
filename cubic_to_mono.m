function [nI, aI, shearI, nII, aII, shearII,...
    bplusI, mplusI, bminusI, mminusI, shapestrainI,...
    bplusII, mplusII, bminusII, mminusII, shapestrainII, Ui, Uj, lambdaI, lambdaII, Q] = cubic_to_mono(a0, a, b, c, phi, CV1, CV2)


%% Transformation matrices
gamma1=(a*(sqrt(2)*a+c*sin(phi)))/(a0*sqrt(2*a^2+c^2+2*sqrt(2)*a*c*sin(phi)));
eps1=a*c*cos(phi)/(sqrt(2)*a0*sqrt(2*a^2+c^2+2*sqrt(2)*a*c*sin(phi)));
alpha1=1/(2*sqrt(2)*a0)*((c*(c+sqrt(2)*a*sin(phi)))/sqrt(2*a^2+c^2+2*sqrt(2)*a*c*sin(phi))+b);
delta1=1/(2*sqrt(2)*a0)*((c*(c+sqrt(2)*a*sin(phi)))/sqrt(2*a^2+c^2+2*sqrt(2)*a*c*sin(phi))-b);

U  = zeros(3,3,12);
U(:,:,1)  = [gamma1 eps1 eps1; eps1 alpha1 delta1; eps1 delta1 alpha1];
U(:,:,2)  = [gamma1 -1*eps1 -1*eps1; -1*eps1 alpha1 delta1; -1*eps1 delta1 alpha1];
U(:,:,3)  = [gamma1 -1*eps1 eps1; -1*eps1 alpha1 -1*delta1; eps1 -1*delta1 alpha1];
U(:,:,4)  = [gamma1 eps1 -1*eps1; eps1 alpha1 -1*delta1; -1*eps1 -1*delta1 alpha1];
U(:,:,5)  = [alpha1 eps1 delta1; eps1 gamma1 eps1; delta1 eps1 alpha1];
U(:,:,6)  = [alpha1 -1*eps1 delta1; -1*eps1 gamma1 -1*eps1; delta1 -1*eps1 alpha1];
U(:,:,7)  = [alpha1 -1*eps1 -1*delta1; -1*eps1 gamma1 eps1; -1*delta1 eps1 alpha1];
U(:,:,8)  = [alpha1 eps1 -1*delta1; eps1 gamma1 -1*eps1; -1*delta1 -1*eps1 alpha1];
U(:,:,9)  = [alpha1 delta1 eps1; delta1 alpha1 eps1; eps1 eps1 gamma1];
U(:,:,10) = [alpha1 delta1 -1*eps1; delta1 alpha1 -1*eps1; -1*eps1 -1*eps1 gamma1];
U(:,:,11) = [alpha1 -1*delta1 eps1; -1*delta1 alpha1 -1*eps1; eps1 -1*eps1 gamma1];
U(:,:,12) = [alpha1 -1*delta1 -1*eps1; -1*delta1 alpha1 eps1; -1*eps1 eps1 gamma1];

ahat=c/(sqrt(2)*a0);  bhat=b/(sqrt(2)*a0);  ghat=a/a0;
Tt2m = [ghat 1/sqrt(2)*ahat*cos(phi) 1/sqrt(2)*ahat*cos(phi);0 0.5*(ahat*sin(phi)+bhat) 0.5*(ahat*sin(phi)-bhat);0 0.5*(ahat*sin(phi)-bhat) 0.5*(ahat*sin(phi)+bhat)];
Q = Tt2m / U(:,:,1);
f1m=a*[1;0;0];  f2m=b/sqrt(2)*[0;-1;1];  f3m=-c/sqrt(2)*[sqrt(2)*cos(phi);sin(phi);sin(phi)]; 
f1m=Q'*f1m;  f2m=Q'*f2m;  f3m=Q'*f3m; 


%% Calculations for variant pair
Ui = U(:,:,CV1);
Uj = U(:,:,CV2);


%%
% First type of twins
kappa=1;
% [nI,aI,etaI,KI,shearI]=twinning_elements_B19p(Ui,Uj,kappa,f1m,f2m,f3m)
[nI,aI,etaI,KI,shearI]=twinning_elements_HS(Ui,Uj,kappa);
[bplusI,mplusI,bminusI,mminusI,lambdaI,shapestrainI]=habit_elements_HS(Uj,nI,aI,1);

% Second type of twins
kappa=-1;
% [nII,aII,etaII,KII,shearII]=twinning_elements_B19p(Ui,Uj,kappa,f1m,f2m,f3m)
[nII,aII,etaII,KII,shearII]=twinning_elements_HS(Ui,Uj,kappa);
[bplusII,mplusII,bminusII,mminusII,lambdaII,shapestrainII]=habit_elements_HS(Uj,nII,aII,2);