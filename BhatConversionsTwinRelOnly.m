% Code to compute the axial strain IPFs of different twin systems
% Written by Aaron Stebner, Northwestern University, 2012
% Modified by Ashley Bucsek, Colorado School of Mines, 2015
% Modified again by ^ in 2016 because want:
% Input = cubicOrientationQuaternion, B2 and B19' lattice parameters
%
% Output = array(192, 18);
% 1:6     CV1    CV2    quatCV1_1 quatCV1_2 quatCV1_3 quatCV1_4   ...
% 7:13    quatCV2_1 quatCV2_2 quatCV2_3 quatCV2_4     m1 m2 m3    ...
% 14:19   n1 n2 n3     LambdaCV1     strain    
%
% Translation of Output ^:
% CV1 is CV i, CV2 j, quat is the quaternion of the monoclinic orientation
% of the CV (first i's then j's), m is the habit plane, n is the twin
% plane, lambda is the twin phase fraction of CV i, strain is the twinning  
% strain in the axial direction (y) from CV j twinning to CV i
% 
% NOTE: m and n are B2 hkls


function Output = BhatConversionsTwinRelOnly(cubicQuat,a0lat,alat,blat,clat,betadeg)

%% Defined Constants
% Lattice Parameters - NiTi
a0 = a0lat;
lat(1) = alat;
lat(2) = blat;
lat(3) = clat;
beta = betadeg*pi/180;  % **********


%% Calculations
CVPairs = [1 2; 3 4; 5 6; 7 8; 9 10; 11 12; ...
    1 3; 1 4; 2 3; 2 4; 5 7; 5 8; 6 7; 6 8; 9 11; 9 12; 10 11; 10 12; ...  % Set 2 / Mode B
    1 5; 1 9; 2 7; 2 12; 3 6; 3 11; 4 8; 4 10; 5 9; 6 11; 7 12; 8 10; ...  % Set 3 / Mode C
    1 8; 1 11; 2 6; 2 10; 3 7; 3 9; 4 5; 4 12; 5 12; 6 10; 7 9; 8 11];

Output = [];
for i = 1 : size(CVPairs, 1)
    [nI, aI, shearI, nII, aII, shearII,...
        bplusI, mplusI, bminusI, mminusI, shapestrainI,...
        bplusII, mplusII, bminusII, mminusII, shapestrainII, Ui, Uj, lamI, lamII] ...
        = cubic_to_mono(a0, lat(1), lat(2), lat(3), beta, CVPairs(i,1), CVPairs(i,2));
    
    [RtwinI, RtwinII] = twinplanerot(Ui, Uj, aI, aII, nI, nII);  % applied on Ui
    
    Rhat = eye(3);
    [quatCV1plusI, quatCV2plusI, R1, U1, R2, U2] = CV2Hexrd(RtwinI, Rhat, Ui, Uj, CVPairs(i,1), CVPairs(i,2), a0, lat(1), lat(2), lat(3), beta);
    [quatCV1plusII, quatCV2plusII, R1, U1, R2, U2] = CV2Hexrd(RtwinII, Rhat, Ui, Uj, CVPairs(i,1), CVPairs(i,2), a0, lat(1), lat(2), lat(3), beta);
    
    % Strains are from CV2 to CV1
    A = quat2rot(cubicQuat);
    F = eye(3) + shearI * dyad(unit(aI), unit(nI));  F = A * F * A';
    strainI_LD = (F*[0;1;0])'*[0;1;0] - 1;
    F = eye(3) + shearII * dyad(unit(aII), unit(nII));  F = A * F * A';
    strainII_LD = (F*[0;1;0])'*[0;1;0] - 1;

    quatCV1plusI = rot2quat( quat2rot(cubicQuat) * quat2rot(quatCV1plusI) );
    quatCV1plusII = rot2quat( quat2rot(cubicQuat) * quat2rot(quatCV1plusII) );
    quatCV2plusI = rot2quat( quat2rot(cubicQuat) * quat2rot(quatCV2plusI) );
    quatCV2plusII = rot2quat( quat2rot(cubicQuat) * quat2rot(quatCV2plusII) );
    Output1 = [CVPairs(i,1) CVPairs(i,2) quatCV1plusI quatCV2plusI mplusI' nI' lamI strainI_LD];
    Output2 = [CVPairs(i,1) CVPairs(i,2) quatCV1plusII quatCV2plusII mplusII' nII' lamII strainII_LD];
    Output = vertcat(Output,Output1,Output2);
    
    
    temp=CVPairs(i,1);  CVPairs(i,1) = CVPairs(i,2);  CVPairs(i,2) = temp;
    [nI, aI, shearI, nII, aII, shearII,...
        bplusI, mplusI, bminusI, mminusI, shapestrainI,...
        bplusII, mplusII, bminusII, mminusII, shapestrainII, Ui, Uj, lamI, lamII] ...
        = cubic_to_mono(a0, lat(1), lat(2), lat(3), beta, CVPairs(i,1), CVPairs(i,2));
    
    [RtwinI, RtwinII] = twinplanerot(Ui, Uj, aI, aII, nI, nII);  % applied on Ui
    
    [quatCV1plusI, quatCV2plusI, R1, U1, R2, U2] = CV2Hexrd(RtwinI, Rhat, Ui, Uj, CVPairs(i,1), CVPairs(i,2), a0, lat(1), lat(2), lat(3), beta);
    [quatCV1plusII, quatCV2plusII, R1, U1, R2, U2] = CV2Hexrd(RtwinII, Rhat, Ui, Uj, CVPairs(i,1), CVPairs(i,2), a0, lat(1), lat(2), lat(3), beta);
    
    A = quat2rot(cubicQuat);
    F = eye(3) + shearI * dyad(unit(aI), unit(nI));  F = A * F * A';
    strainI_LD = (F*[0;1;0])'*[0;1;0] - 1;
    F = eye(3) + shearII * dyad(unit(aII), unit(nII));  F = A * F * A';
    strainII_LD = (F*[0;1;0])'*[0;1;0] - 1;    

    quatCV1plusI = rot2quat( quat2rot(cubicQuat) * quat2rot(quatCV1plusI) );
    quatCV1plusII = rot2quat( quat2rot(cubicQuat) * quat2rot(quatCV1plusII) );
    quatCV2plusI = rot2quat( quat2rot(cubicQuat) * quat2rot(quatCV2plusI) );
    quatCV2plusII = rot2quat( quat2rot(cubicQuat) * quat2rot(quatCV2plusII) );
    Output1 = [CVPairs(i,1) CVPairs(i,2) quatCV1plusI quatCV2plusI mplusI' nI' lamI strainI_LD];
    Output2 = [CVPairs(i,1) CVPairs(i,2) quatCV1plusII quatCV2plusII mplusII' nII' lamII strainII_LD];
    Output = vertcat(Output,Output1,Output2);
end