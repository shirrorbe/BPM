function [w] = SlerpBlendMobius(z ,zi, zj, zk, a_ijk,...
    b_ijk,c_ijk,d_ijk, a_ij,b_ij,c_ij,d_ij, a_jk,b_jk,c_jk,d_jk, a_ki,b_ki,c_ki,d_ki)
% function [w] = SlerpBlendMobius(z ,zi, zj, zk, a_ijk,...
%     b_ijk,c_ijk,d_ijk, a_ij,b_ij,c_ij,d_ij, a_jk,b_jk,c_jk,d_jk, a_ki,b_ki,c_ki,d_ki)

%BarycentricMobiusBlend Summary of this function goes here
%   Inputs
% z - the internal point
% zi, zj, zk - vertex locations of triangle ijk
% Mij, Mjk, Mki - Mobius of faces near edges eij,ejk,eki

% TODO: make more efficient (dont compute inverse for each point in z, in
% addition can use closed form of inversion of 2X2 matrix
% if inpolygon(real(z),imag(z),real([zi,zj,zk]),imag([zi,zj,zk]))&&real(z)>0.2&&imag(z)>0.9
%     bla=1;
% end

Mijk = [a_ijk b_ijk; c_ijk d_ijk];

if isnan(a_ij) || isnan(b_ij) || isnan(c_ij) || isnan(d_ij)
    Mij = Mijk;
else
    Mij = [a_ij b_ij; c_ij d_ij];
end

if isnan(a_jk) || isnan(b_jk) || isnan(c_jk) || isnan(d_jk)
    Mjk = Mijk;
else
    Mjk = [a_jk b_jk; c_jk d_jk];
end

if isnan(a_ki) || isnan(b_ki) || isnan(c_ki) || isnan(d_ki)
    Mki = Mijk;
else
    Mki = [a_ki b_ki; c_ki d_ki];
end

% Mij = [a_ij b_ij; c_ij d_ij];
% Mjk = [a_jk b_jk; c_jk d_jk];
% Mki = [a_ki b_ki; c_ki d_ki];



r_ij = PointToEdgeDistance(z,zi,zj);
r_jk = PointToEdgeDistance(z,zj,zk);
r_ki = PointToEdgeDistance(z,zk,zi);

eps = 1e-20;
% if r_ij < eps 
%     gamma_ij = 1; gamma_jk = 0; gamma_ki = 0;
% elseif r_jk < eps
%     gamma_ij = 0; gamma_jk = 1; gamma_ki = 0;
% elseif r_ki < eps
%     gamma_ij = 0; gamma_jk = 0; gamma_ki = 1;
% else
%     gamma_ij = 1/r_ij;
%     gamma_jk = 1/r_jk;
%     gamma_ki = 1/r_ki;
% end




if abs(z-zi)<eps
    gamma_ij = 1; gamma_jk = 0; gamma_ki = 0;
elseif abs(z-zj)<eps
    gamma_ij = 0; gamma_jk = 1; gamma_ki = 0;
elseif abs(z-zk)<eps
    gamma_ij = 0; gamma_jk = 0; gamma_ki = 1;
else
    gamma_ij = (r_jk*r_ki)/(r_jk*r_ki+r_ij*r_jk+r_ij*r_ki);
    gamma_jk = (r_ij*r_ki)/(r_jk*r_ki+r_ij*r_jk+r_ij*r_ki);
    gamma_ki = (r_ij*r_jk)/(r_jk*r_ki+r_ij*r_jk+r_ij*r_ki);
end





Mijk_inv = Mijk^-1;

% compute the error term E
if ~isnan(Mij) 
    Eij_p = sqrtm(Mij*Mijk_inv);
    Eij_m = sqrtm(-Mij*Mijk_inv);
    norm_p = norm(Eij_p-eye(2));
    norm_m = norm(Eij_m-eye(2));
    Eij = Eij_p;
    if norm_p > norm_m
        Eij = Eij_m;
    end
else
    Eij = eye(2);
    gamma_ij = 0;
end

if ~isnan(Mjk) 
    Ejk_p = sqrtm(Mjk*Mijk_inv);
    Ejk_m = sqrtm(-Mjk*Mijk_inv);
    norm_p = norm(Ejk_p-eye(2));
    norm_m = norm(Ejk_m-eye(2));
    Ejk = Ejk_p;
    if norm_p > norm_m
        Ejk = Ejk_m;
    end
else
     Ejk = eye(2);
     gamma_jk = 0;
end

if ~isnan(Mki) 
    Eki_p = sqrtm(Mki*Mijk_inv);
    Eki_m = sqrtm(-Mki*Mijk_inv);
    norm_p = norm(Eki_p-eye(2));
    norm_m = norm(Eki_m-eye(2));
    Eki = Eki_p;
    if norm_p > norm_m
        Eki = Eki_m;
    end
else
    Eki = eye(2);
    gamma_ki = 0;
end


% % edge lengths weights
% l_ij = abs(zj-zi); gamma_ij = gamma_ij*l_ij;
% l_jk = abs(zk-zj); gamma_jk = gamma_jk*l_jk;
% l_ki = abs(zk-zi); gamma_ki = gamma_ki*l_ki;
% 


Eijk = gamma_ij*logm(Eij)+gamma_jk*logm(Ejk)+gamma_ki*logm(Eki);

% Eijk = (gamma_ij*logm(Eij)+gamma_jk*logm(Ejk)+gamma_ki*logm(Eki))/(gamma_ij+gamma_jk+gamma_ki);

M = expm(Eijk)*Mijk;

% a_blend = M(1,1); b_blend = M(1,2); c_blend = M(2,1); d_blend = M(2,2);
w = MobiusTransform(z,M(1,1),M(1,2),M(2,1),M(2,2));
end


function d = PointToEdgeDistance(z, z1, z2)
      pt = [real(z),imag(z),0];       
      v1 = [real(z1),imag(z1),0];
      v2 = [real(z2),imag(z2),0];
      a = v1 - v2;
      b = pt - v2;
      d = norm(cross(a,b)) / norm(a);
end