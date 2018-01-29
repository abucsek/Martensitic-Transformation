% Created by Ashley Bucsek, Colorado School of Mines, 2014
% Returns all the hkls in a family of hkls, assuming cubic symmetry

function [vecout,L]=cubic_symmetries(vecin)

e1=[1;0;0]; e2=[0;1;0]; e3=[0;0;1];

rots=cell(24,1);
rots{1,1}=eye(3);
rots{2,1}=Rodrot(pi/2,e1);
rots{3,1}=Rodrot(pi,e1);
rots{4,1}=Rodrot(270*pi/180,e1);
rots{5,1}=Rodrot(pi/2,e2);
rots{6,1}=Rodrot(pi,e2);
rots{7,1}=Rodrot(270*pi/180,e2);
rots{8,1}=Rodrot(pi/2,e3);
rots{9,1}=Rodrot(pi,e3);
rots{10,1}=Rodrot(270*pi/180,e3);
rots{11,1}=Rodrot(pi,unit(e1+e2));
rots{12,1}=Rodrot(pi,unit(e1-e2));
rots{13,1}=Rodrot(pi,unit(e2+e3));
rots{14,1}=Rodrot(pi,unit(e2-e3));
rots{15,1}=Rodrot(pi,unit(e1+e3));
rots{16,1}=Rodrot(pi,unit(e1-e3));
rots{17,1}=Rodrot(120*pi/180,unit(e1+e2+e3));
rots{18,1}=Rodrot(120*pi/180,unit(e1+e2-e3));
rots{19,1}=Rodrot(120*pi/180,unit(e1-e2+e3));
rots{20,1}=Rodrot(120*pi/180,unit(-e1+e2+e3));
rots{21,1}=Rodrot(240*pi/180,unit(e1+e2+e3));
rots{22,1}=Rodrot(240*pi/180,unit(e1+e2-e3));
rots{23,1}=Rodrot(240*pi/180,unit(e1-e2+e3));
rots{24,1}=Rodrot(240*pi/180,unit(-e1+e2+e3));

vecout=zeros(24,3);
for i=1:24
    vecout(i,:)=transpose(rots{i,1}*vecin);
end

i=1;
while i<=size(vecout,1)
    for j=1:3
        if abs(vecout(i,j))<=1e-10
            vecout(i,j)=0;
        end
    end
    i=i+1;
end

vecout=unique(vecout,'rows');

i=1;
while i<=size(vecout,1)
    for j=1:3
        vecout(i,j)=round(vecout(i,j)*1e4)/1e4;
    end
    i=i+1;
end

vecout=unique(vecout,'rows');

L=size(vecout,1);