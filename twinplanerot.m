function [RtwinI, RtwinII] = twinplanerot(Ui, Uj, aI, aII, nI, nII)
% Based on H&S '99 terminology:
% Rij*Ui - Uj = dyad(a, nhat)

RtwinI = (dyad(aI,nI) + Uj) * inv(Ui);
RtwinII = (dyad(aII,nII) + Uj) * inv(Ui);