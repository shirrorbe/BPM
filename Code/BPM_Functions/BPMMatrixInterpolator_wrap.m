function [w,t_weights_ijkt,Mijk,Mij,Mjk,Mki] = BPMMatrixInterpolator_wrap(x,Vz, Vw,f)
% function [w,t_weights_ijkt] = BPMMatrixInterpolator_wrap(x,Vz, Vw,f)

%UNTITLED Summary of this function goes here
% Vz = [vi;vj;vk;vl;vm;vn];
% F = [i j k; k j l; m i k; i n j]; (F = [1 2 3; 3 2 4; 5 1 3; 1 6 2])

switch f
    case 2
%         disp('2')
        Vz =  [Vz(3,:); Vz(2,:); Vz(4,:); nan(2,3); Vz(1,:)];
        Vw =  [Vw(3,:); Vw(2,:); Vw(4,:); nan(2,3); Vw(1,:)];

    case 3
%         disp('3')
        Vz =  [Vz(5,:); Vz(1,:); Vz(3,:); nan(2,3); Vz(2,:)];
        Vw =  [Vw(5,:); Vw(1,:); Vw(3,:); nan(2,3); Vw(2,:)];
    case 4
%        disp('4')
        Vz =  [Vz(1,:); Vz(6,:); Vz(2,:); nan(1,3); Vz(3,:);  nan(1,3)];
        Vw =  [Vw(1,:); Vw(6,:); Vw(2,:); nan(1,3); Vw(3,:); nan(1,3)];

end


zi = complex(Vz(1,1),Vz(1,2)); zj = complex(Vz(2,1),Vz(2,2));
zk = complex(Vz(3,1),Vz(3,2)); zl = complex(Vz(4,1),Vz(4,2));
zm = complex(Vz(5,1),Vz(5,2)); zn = complex(Vz(6,1),Vz(6,2));

zf_ijk = [zi, zj, zk]; zf_ij = [zi, zn, zj]; zf_jk = [zk, zj, zl];
zf_ki= [zm, zi, zk];

wi = complex(Vw(1,1),Vw(1,2)); wj = complex(Vw(2,1),Vw(2,2));
wk = complex(Vw(3,1),Vw(3,2)); wl = complex(Vw(4,1),Vw(4,2));
wm = complex(Vw(5,1),Vw(5,2)); wn = complex(Vw(6,1),Vw(6,2));

wf_ijk = [wi, wj, wk]; wf_ij = [wi, wn, wj]; wf_jk = [wk, wj, wl];
wf_ki= [wm, wi, wk];


% 
% zf_ijk = [complex(Vz(1,1),Vz(1,2)), complex(Vz(2,1),Vz(2,2)), complex(Vz(3,1),Vz(3,2))];
% zf_ij = [complex(Vz(1,1),Vz(1,2)), complex(Vz(6,1),Vz(6,2)), complex(Vz(2,1),Vz(2,2))];
% zf_jk = [complex(Vz(3,1),Vz(3,2)), complex(Vz(2,1),Vz(2,2)), complex(Vz(4,1),Vz(4,2))];
% zf_ki = [complex(Vz(5,1),Vz(5,2)), complex(Vz(1,1),Vz(1,2)), complex(Vz(3,1),Vz(3,2))];
% 
% wf_ijk = [complex(Vw(1,1),Vw(1,2)), complex(Vw(2,1),Vw(2,2)), complex(Vw(3,1),Vw(3,2))];
% wf_ij = [complex(Vw(1,1),Vw(1,2)), complex(Vw(6,1),Vw(6,2)), complex(Vw(2,1),Vw(2,2))];
% wf_jk = [complex(Vw(3,1),Vw(3,2)), complex(Vw(2,1),Vw(2,2)), complex(Vw(4,1),Vw(4,2))];
% wf_ki = [complex(Vw(5,1),Vw(5,2)), complex(Vw(1,1),Vw(1,2)), complex(Vw(3,1),Vw(3,2))];



z = complex(x(:,1),x(:,2));

[a_ijk,b_ijk,c_ijk,d_ijk] = ComputeMobiusCoeffsFromPos(zf_ijk(1), zf_ijk(2), zf_ijk(3),  wf_ijk(1), wf_ijk(2), wf_ijk(3));  

[a_ij,b_ij,c_ij,d_ij] = arrayfun(@ComputeMobiusCoeffsFromPos,  zf_ij(1), zf_ij(2), zf_ij(3), wf_ij(1), wf_ij(2), wf_ij(3));  

[a_jk,b_jk,c_jk,d_jk] = arrayfun(@ComputeMobiusCoeffsFromPos,  zf_jk(1), zf_jk(2), zf_jk(3), wf_jk(1), wf_jk(2), wf_jk(3));  

[a_ki,b_ki,c_ki,d_ki] = arrayfun(@ComputeMobiusCoeffsFromPos,  zf_ki(1), zf_ki(2), zf_ki(3), wf_ki(1), wf_ki(2), wf_ki(3));  

Mijk = [a_ijk b_ijk; c_ijk d_ijk];
Mij = [a_ij b_ij; c_ij d_ij];
Mjk = [a_jk b_jk; c_jk d_jk];
Mki = [a_ki b_ki; c_ki d_ki];
[O_mats,weights_ijkt,Mijk,Mij,Mjk,Mki] = BPMMatrixInterpolator(z,zf_ijk,Mijk,Mij,Mjk,Mki);

% [O_mats,weights_ijkt] = BPMMatrixInterpolator(z,zf_ijk,Mijk,Mij,Mjk,Mki);
w = MobiusTransform(z,O_mats(:,1),O_mats(:,2),O_mats(:,3),O_mats(:,4));
% t_weights_ijkt = weights_ijkt(:,3); 
t_weights_ijkt = weights_ijkt; 

end