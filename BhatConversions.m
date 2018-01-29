% Calculate all 192 HPV pairs using Crystallographic Theory of Martensite
% eigenvalue calculations, convert these to orientations relative to
% standard a||X, c*||Z B19' orientations
%
% Written by Ashley Bucsek, Colorado School of Mines, 2016
%
% Input = cubicOrientationQuaternion, B2 and B19' lattice parameters
%
% Output = array(192, 18);
% 1:6     CV1    CV2    quatCV1_1 quatCV1_2 quatCV1_3 quatCV1_4   ...
% 7:13    quatCV2_1 quatCV2_2 quatCV2_3 quatCV2_4     m1 m2 m3    ...
% 14:19   n1 n2 n3     LambdaCV1     strain    stressStrainWork
%
% Translation of Output ^:
% CV1 is CV i, CV2 j, quat is the quaternion of the monoclinic orientation
% of the CV (first i's then j's), m is the habit plane, n is the twin
% plane, lambda is the twin phase fraction of CV i, strain is the transformation 
% strain in the axial direction (y), and stressStrainWork is the work
% produced
% 
% NOTE: m and n are B2 hkls


function [Output, Bhat] = BhatConversions(cubicQuat,a0lat,alat,blat,clat,betadeg)


%% Defined Constants
% Lattice Parameters - NiTi
a0 = a0lat;
lat(1) = alat;
lat(2) = blat;
lat(3) = clat;
beta = betadeg*pi/180;

% Assume a uniaxial load in the y-direction (normalized)
% (Update with actuall stress tensor if known)
Stress = [0 0 0; 0 1 0; 0 0 0];


%% Calculations
% Considering compatible habit plane twins only
CVPairs = [1 3; 1 4; 2 3; 2 4; 5 7; 5 8; 6 7; 6 8; 9 11; 9 12; 10 11; 10 12; ...  % Set 2 / Mode B
    1 5; 1 9; 2 7; 2 12; 3 6; 3 11; 4 8; 4 10; 5 9; 6 11; 7 12; 8 10];            % Set 3 / Mode C

Output = [];
for i = 1 : size(CVPairs, 1)

    [nI, aI, shearI, nII, aII, shearII,...
        bplusI, mplusI, bminusI, mminusI, shapestrainI,...
        bplusII, mplusII, bminusII, mminusII, shapestrainII, Ui, Uj, lamI, lamII, Q] ...
        = cubic_to_mono(a0, lat(1), lat(2), lat(3), beta, CVPairs(i,1), CVPairs(i,2));
    
    [RtwinI, RtwinII] = twinplanerot(Ui, Uj, aI, aII, nI, nII);  % applied on Ui
    [RhatplusI, RhatminusI, RhatplusII, RhatminusII] = habitplanerots(Uj, lamI, lamII, aI, aII, nI, nII, ...
        bplusI, bminusI, bplusII, bminusII, mplusI, mminusI, mplusII, mminusII);  % applied on (Uj + lam*dyad(a,n))
    
    
    [quatCV1plusI, quatCV2plusI, R1, U1, R2, U2] = CV2Hexrd(RtwinI, RhatplusI, Ui, Uj, CVPairs(i,1), CVPairs(i,2), a0, lat(1), lat(2), lat(3), beta, Q);
    [quatCV1minusI, quatCV2minusI, R1, U1, R2, U2] = CV2Hexrd(RtwinI, RhatminusI, Ui, Uj, CVPairs(i,1), CVPairs(i,2), a0, lat(1), lat(2), lat(3), beta, Q);
    [quatCV1plusII, quatCV2plusII, R1, U1, R2, U2] = CV2Hexrd(RtwinII, RhatplusII, Ui, Uj, CVPairs(i,1), CVPairs(i,2), a0, lat(1), lat(2), lat(3), beta, Q);
    [quatCV1minusII, quatCV2minusII, R1, U1, R2, U2] = CV2Hexrd(RtwinII, RhatminusII, Ui, Uj, CVPairs(i,1), CVPairs(i,2), a0, lat(1), lat(2), lat(3), beta, Q);
    
    A = quat2rot(cubicQuat);
    F = eye(3) + shapestrainI * dyad(unit(bplusI), unit(mplusI));  F = A * F * A';
    strainplusI_LD = (1/2*(F'*F-eye(3))*[0;1;0])'*[0;1;0] ;
    stressStrainWork_plusI = doubleDot(Stress,1/2*(F'*F-eye(3)));
    F = eye(3) + shapestrainI * dyad(unit(bminusI), unit(mminusI));  F = A * F * A';
    strainminusI_LD = (1/2*(F'*F-eye(3))*[0;1;0])'*[0;1;0] ;
    stressStrainWork_minusI = doubleDot(Stress,1/2*(F'*F-eye(3)));
    F = eye(3) + shapestrainII * dyad(unit(bplusII), unit(mplusII));  F = A * F * A';
    strainplusII_LD = (1/2*(F'*F-eye(3))*[0;1;0])'*[0;1;0] ;
    stressStrainWork_plusII = doubleDot(Stress,1/2*(F'*F-eye(3)));
    F = eye(3) + shapestrainII * dyad(unit(bminusII), unit(mminusII));  F = A * F * A';
    strainminusII_LD = (1/2*(F'*F-eye(3))*[0;1;0])'*[0;1;0] ;
    stressStrainWork_minusII = doubleDot(Stress,1/2*(F'*F-eye(3)));
    
    quatCV1plusI = rot2quat( quat2rot(cubicQuat) * quat2rot(quatCV1plusI) );
    quatCV1minusI = rot2quat( quat2rot(cubicQuat) * quat2rot(quatCV1minusI) );
    quatCV1plusII = rot2quat( quat2rot(cubicQuat) * quat2rot(quatCV1plusII) );
    quatCV1minusII = rot2quat( quat2rot(cubicQuat) * quat2rot(quatCV1minusII) );
    quatCV2plusI = rot2quat( quat2rot(cubicQuat) * quat2rot(quatCV2plusI) );
    quatCV2minusI = rot2quat( quat2rot(cubicQuat) * quat2rot(quatCV2minusI) );
    quatCV2plusII = rot2quat( quat2rot(cubicQuat) * quat2rot(quatCV2plusII) );
    quatCV2minusII = rot2quat( quat2rot(cubicQuat) * quat2rot(quatCV2minusII) );    
    Output1 = [CVPairs(i,1) CVPairs(i,2) quatCV1plusI quatCV2plusI mplusI' nI' lamI strainplusI_LD stressStrainWork_plusI];
    Output2 = [CVPairs(i,1) CVPairs(i,2) quatCV1minusI quatCV2minusI mminusI' nI' lamI strainminusI_LD stressStrainWork_minusI];
    Output3 = [CVPairs(i,1) CVPairs(i,2) quatCV1plusII quatCV2plusII mplusII' nII' lamII strainplusII_LD stressStrainWork_plusII];
    Output4= [CVPairs(i,1) CVPairs(i,2) quatCV1minusII quatCV2minusII mminusII' nII' lamII strainminusII_LD stressStrainWork_minusII];
    Output = vertcat(Output,Output1,Output2,Output3,Output4);

    
    % Switch CV1 & CV2 to find other half of solutions
    temp=CVPairs(i,1);  CVPairs(i,1) = CVPairs(i,2);  CVPairs(i,2) = temp;
    
    [nI, aI, shearI, nII, aII, shearII,...
        bplusI, mplusI, bminusI, mminusI, shapestrainI,...
        bplusII, mplusII, bminusII, mminusII, shapestrainII, Ui, Uj, lamI, lamII, Q] ...
        = cubic_to_mono(a0, lat(1), lat(2), lat(3), beta, CVPairs(i,1), CVPairs(i,2));
    
    [RtwinI, RtwinII] = twinplanerot(Ui, Uj, aI, aII, nI, nII);  % applied on Ui
    [RhatplusI, RhatminusI, RhatplusII, RhatminusII] = habitplanerots(Uj, lamI, lamII, aI, aII, nI, nII, ...
        bplusI, bminusI, bplusII, bminusII, mplusI, mminusI, mplusII, mminusII);  % applied on (Uj + lam*dyad(a,n))
    
    [quatCV1plusI, quatCV2plusI, R1, U1, R2, U2] = CV2Hexrd(RtwinI, RhatplusI, Ui, Uj, CVPairs(i,1), CVPairs(i,2), a0, lat(1), lat(2), lat(3), beta, Q);
    [quatCV1minusI, quatCV2minusI, R1, U1, R2, U2] = CV2Hexrd(RtwinI, RhatminusI, Ui, Uj, CVPairs(i,1), CVPairs(i,2), a0, lat(1), lat(2), lat(3), beta, Q);
    [quatCV1plusII, quatCV2plusII, R1, U1, R2, U2] = CV2Hexrd(RtwinII, RhatplusII, Ui, Uj, CVPairs(i,1), CVPairs(i,2), a0, lat(1), lat(2), lat(3), beta, Q);
    [quatCV1minusII, quatCV2minusII, R1, U1, R2, U2] = CV2Hexrd(RtwinII, RhatminusII, Ui, Uj, CVPairs(i,1), CVPairs(i,2), a0, lat(1), lat(2), lat(3), beta, Q);
    
    A = quat2rot(cubicQuat);
    F = eye(3) + shapestrainI * dyad(unit(bplusI), unit(mplusI));  F = A * F * A';
    strainplusI_LD = (1/2*(F'*F-eye(3))*[0;1;0])'*[0;1;0] ;
    stressStrainWork_plusI = doubleDot(Stress,1/2*(F'*F-eye(3)));
    F = eye(3) + shapestrainI * dyad(unit(bminusI), unit(mminusI));  F = A * F * A';
    strainminusI_LD = (1/2*(F'*F-eye(3))*[0;1;0])'*[0;1;0] ;
    stressStrainWork_minusI = doubleDot(Stress,1/2*(F'*F-eye(3)));
    F = eye(3) + shapestrainII * dyad(unit(bplusII), unit(mplusII));  F = A * F * A';
    strainplusII_LD = (1/2*(F'*F-eye(3))*[0;1;0])'*[0;1;0] ;
    stressStrainWork_plusII = doubleDot(Stress,1/2*(F'*F-eye(3)));
    F = eye(3) + shapestrainII * dyad(unit(bminusII), unit(mminusII));  F = A * F * A';
    strainminusII_LD = (1/2*(F'*F-eye(3))*[0;1;0])'*[0;1;0] ;
    stressStrainWork_minusII = doubleDot(Stress,1/2*(F'*F-eye(3)));
    
    quatCV1plusI = rot2quat( quat2rot(cubicQuat) * quat2rot(quatCV1plusI) );
    quatCV1minusI = rot2quat( quat2rot(cubicQuat) * quat2rot(quatCV1minusI) );
    quatCV1plusII = rot2quat( quat2rot(cubicQuat) * quat2rot(quatCV1plusII) );
    quatCV1minusII = rot2quat( quat2rot(cubicQuat) * quat2rot(quatCV1minusII) );
    quatCV2plusI = rot2quat( quat2rot(cubicQuat) * quat2rot(quatCV2plusI) );
    quatCV2minusI = rot2quat( quat2rot(cubicQuat) * quat2rot(quatCV2minusI) );
    quatCV2plusII = rot2quat( quat2rot(cubicQuat) * quat2rot(quatCV2plusII) );
    quatCV2minusII = rot2quat( quat2rot(cubicQuat) * quat2rot(quatCV2minusII) );    
    Output1 = [CVPairs(i,1) CVPairs(i,2) quatCV1plusI quatCV2plusI mplusI' nI' lamI strainplusI_LD stressStrainWork_plusI];
    Output2 = [CVPairs(i,1) CVPairs(i,2) quatCV1minusI quatCV2minusI mminusI' nI' lamI strainminusI_LD stressStrainWork_minusI];
    Output3 = [CVPairs(i,1) CVPairs(i,2) quatCV1plusII quatCV2plusII mplusII' nII' lamII strainplusII_LD stressStrainWork_plusII];
    Output4= [CVPairs(i,1) CVPairs(i,2) quatCV1minusII quatCV2minusII mminusII' nII' lamII strainminusII_LD stressStrainWork_minusII];
    Output = vertcat(Output,Output1,Output2,Output3,Output4);
    
end

[temp, ind] = sort(Output(:,18), 'descend');
Output = Output(ind,:);