function [w] = SlerpBlendMobius(z ,zi, zj, zk, a_ijk,...
    b_ijk,c_ijk,d_ijk, a_ij,b_ij,c_ij,d_ij, a_jk,b_jk,c_jk,d_jk, a_ki,b_ki,c_ki,d_ki)

%BarycentricMobiusBlend Summary of this function goes here
%   Inputs
% z - the internal point
% zi, zj, zk - vertex locations of triangle ijk
% Mij, Mjk, Mki - Mobius of faces near edges eij,ejk,eki


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


r_ij = PointToEdgeDistance(z,zi,zj);
r_jk = PointToEdgeDistance(z,zj,zk);
r_ki = PointToEdgeDistance(z,zk,zi);

eps = 1e-20;

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
    delta_ij = Mij*Mijk_inv;
    Eij = sign(trace(real(delta_ij)))*delta_ij;
else
    Eij = eye(2);
    gamma_ij = 0;
end

if ~isnan(Mjk) 
    delta_jk = Mjk*Mijk_inv;
    Ejk = sign(trace(real(delta_jk)))*delta_jk;
else
    Ejk = eye(2);
    gamma_jk = 0;
end

if ~isnan(Mki) 
    delta_ki =Mki*Mijk_inv; 
    Eki = sign(trace(real(delta_ki)))*delta_ki;

else
    Eki = eye(2);
    gamma_ki = 0;
end




w_t = 0.5;
Eijk = (1-w_t)*gamma_ij*logm(Eij)+(1-w_t)*gamma_jk*logm(Ejk)+(1-w_t)*gamma_ki*logm(Eki);



M = expm(Eijk)*Mijk;


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
