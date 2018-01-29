function [quatCV1, quatCV2, R1, U1, R2, U2] = CV2Hexrd(R, Rhat, Ui, Uj, CVi, CVj, a0, a, b, c, beta, Q)
e1 = [1;0;0];
e2 = [0;1;0];
e3 = [0;0;1];

c1 = a0 * [1;0;0];
c2 = a0 * [0;1;0];
c3 = a0 * [0;0;1];

for ii = 1:2
    
    if ii == 1
        CV = CVi;   U_CV = Ui;
    elseif ii == 2
        CV = CVj;   U_CV = Uj;
    end
    
    %% Rotated (by B2 orientation) tetragonal
    if CV < 5
        % CV's 1-4
        TetragonalStretch = [1 0 0; 0 sqrt(2) 0; 0 0 sqrt(2)];
        R45 = Rodrot(-45*pi/180, e1);
        Tetragonal = R45 * TetragonalStretch;
    elseif CV < 9
        % CV's 5-8
        TetragonalStretch = [sqrt(2) 0 0; 0 1 0; 0 0 sqrt(2)];
        R45 = Rodrot(-45*pi/180, e2);
        Tetragonal = R45 * TetragonalStretch;
    else
        % CV's 9-12
        TetragonalStretch = [sqrt(2) 0 0; 0 sqrt(2) 0; 0 0 1];
        R45 = Rodrot(-45*pi/180, e3);
        Tetragonal = R45 * TetragonalStretch;
    end
    
    t1 = Tetragonal*c1;  t2 = Tetragonal*c2;  t3 = Tetragonal*c3;
    
    %% Rotated (by B2 orientation) monoclinic - initial
        if ii == 1
            m1 = Rhat * R * U_CV * t1;
            m2 = Rhat * R * U_CV * t2;
            m3 = Rhat * R * U_CV * t3;
        elseif ii == 2
            m1 = Rhat * U_CV * t1;
            m2 = Rhat * U_CV * t2;
            m3 = Rhat * U_CV * t3;
        end
    M = [m1 m2 m3];
    
    % FIX #1
    % Fixing order of lattice vectors s.t. mono = [m1 m2 m3]
    CHECKS = [norm(m1) norm(m2) norm(m3);
        ang(m1,m2)*180/pi ang(m2,m3)*180/pi ang(m3,m1)*180/pi];
    [temp ind] = sort(CHECKS(1,:),2);
    m1 = M(:,ind(1));
    m2 = M(:,ind(2));
    m3 = M(:,ind(3));
    
    % FIX #2
    % Fixing direction of m3 s.t. right-handed orthogonal system
    if ang(m3,m1)*180/pi < 90
        m3 = -m3;
    end
    
    %Calculate R via polar decomp.
    M_unrotated = [a 0 c*cos(beta);
        0 b 0;
        0 0 c*sin(beta)];
    M_rotated = [m1 m2 m3];
    F = M_rotated / M_unrotated;
    [RHexrd, Utemp] = polardecomp(F);
    
    % FIX #3
    % Remove inversion from orientation
    if det(RHexrd) < 0
        RHexrd = -RHexrd;
    end
    
    if ii == 1
        quatCV1 = rot2quat(RHexrd);
        R1 = RHexrd;
        U1 = Utemp;
    elseif ii == 2
        quatCV2 = rot2quat(RHexrd);
        R2 = RHexrd;
        U2 = Utemp;
    end
end