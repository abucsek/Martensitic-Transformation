function [RhatplusI, RhatminusI, RhatplusII, RhatminusII] = ...
    habitplanerots(Uj, lamI, lamII, aI, aII, nI, nII, bplusI, bminusI, bplusII, bminusII, mplusI, mminusI, mplusII, mminusII)
% Based on H&S '99 terminology:
% Rhatij * (lambda*Rij*Ui + (1-lambda)*Uj) - I = dyad(b,mhat)
% Rhatij * (Uj + lambda*dyad(a,nhat)) - I = dyad(b,mhat)

RhatplusI = (dyad(bplusI,mplusI) + eye(3)) / (Uj + lamI * dyad(aI,nI));
RhatminusI = (dyad(bminusI,mminusI) + eye(3)) / (Uj + lamI * dyad(aI,nI));

RhatplusII = (dyad(bplusII,mplusII) + eye(3)) / (Uj + lamII * dyad(aII,nII));
RhatminusII = (dyad(bminusII,mminusII) + eye(3)) / (Uj + lamII * dyad(aII,nII));