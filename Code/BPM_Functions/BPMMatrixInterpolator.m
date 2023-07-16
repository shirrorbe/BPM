function [O_mats,weights_ijkt,Mijk,Mij,Mjk,Mki] = BPMMatrixInterpolator(z,z_ijk,Mijk,Mij,Mjk,Mki)
% function [O_mats,weights_ijkt] = BPMMatrixInterpolator(z,z_ijk,Mijk,Mij,Mjk,Mki)

%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% a_ijk = Mijk(1,1);b_ijk = Mijk(1,2);c_ijk = Mijk(2,1);d_ijk = Mijk(2,2);
a_ij = Mij(1,1);b_ij = Mij(1,2);c_ij = Mij(2,1);d_ij = Mij(2,2);
a_jk = Mjk(1,1);b_jk = Mjk(1,2);c_jk = Mjk(2,1);d_jk = Mjk(2,2);
a_ki = Mki(1,1);b_ki = Mki(1,2);c_ki = Mki(2,1);d_ki = Mki(2,2);

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







% compute the error term E
Mijk_inv = Mijk^-1;

if ~isnan(Mij) 
    delta_ij = Mij*Mijk_inv;
    Eij = sign(trace(real(delta_ij)))*delta_ij;
else
    Eij = eye(2);
end

if ~isnan(Mjk) 
    delta_jk = Mjk*Mijk_inv;
    Ejk = sign(trace(real(delta_jk)))*delta_jk;
else
    Ejk = eye(2);
end

if ~isnan(Mki) 
    delta_ki =Mki*Mijk_inv; 
    Eki = sign(trace(real(delta_ki)))*delta_ki;

else
    Eki = eye(2);
end



[gamma_ij,gamma_jk,gamma_ki,r_ij,r_jk,r_ki] = BarycentricWeights(z,z_ijk);

log_Eij = logm(Eij);log_Ejk = logm(Ejk);log_Eki = logm(Eki);
log_Eij=repmat(reshape(log_Eij,1,[]),size(gamma_ij,1),1);
log_Ejk=repmat(reshape(log_Ejk,1,[]),size(gamma_ij,1),1);
log_Eki=repmat(reshape(log_Eki,1,[]),size(gamma_ij,1),1);

gamma_ij_vec = repmat(gamma_ij,1,4);
gamma_jk_vec = repmat(gamma_jk,1,4);
gamma_ki_vec = repmat(gamma_ki,1,4);

% Eijk_vec = gamma_ij_vec.*log_Eij+gamma_jk_vec.*log_Ejk+gamma_ki_vec.*log_Eki;


w_t = 0.5;
Eijk_vec = (1-w_t)*gamma_ij_vec.*log_Eij+(1-w_t)*gamma_jk_vec.*log_Ejk+(1-w_t)*gamma_ki_vec.*log_Eki;
weights_ijkt = [(1-w_t)*gamma_ij,(1-w_t)*gamma_jk,(1-w_t)*gamma_ki,w_t*ones(size(gamma_ij))];

% 
% zi = z_ijk(1);zj = z_ijk(2);zk = z_ijk(3);
% zm = (zi+zj+zk)/3;
% rm = max(PointToEdgeDistance(zm,zi,zj),max(PointToEdgeDistance(zm,zi,zj),PointToEdgeDistance(zm,zi,zj)));
% min_r = min(min(r_ij,r_jk),r_ki);
% w_t = -0.2*min_r./rm + 0.5;
% % w_t = -0.3*min_r./rm + 0.5;
% % w_t = 0.2*min_r./rm + 0.5;
% Eijk_vec = (1-w_t).*gamma_ij_vec.*log_Eij+(1-w_t).*gamma_jk_vec.*log_Ejk+(1-w_t).*gamma_ki_vec.*log_Eki;
% weights_ijkt = [(1-w_t).*gamma_ij,(1-w_t).*gamma_jk,(1-w_t).*gamma_ki,w_t.*ones(size(gamma_ij))];
% 



O_mats = zeros(size(Eijk_vec));
for ii=1:size(gamma_ij,1)
    Eijk = reshape(Eijk_vec(ii,:),2,2);
    O_mat = expm(Eijk)*Mijk;
    O_mats(ii,:) = O_mat(:);
end

O_mats = O_mats(:,[1,3,2,4]);

% O_mat = expm(Eijk)*Mijk;



end


function [gamma_ij,gamma_jk,gamma_ki,r_ij,r_jk,r_ki] = BarycentricWeights(z,z_ijk)

r_ij = PointToEdgeDistance(z,z_ijk(1),z_ijk(2));
r_jk = PointToEdgeDistance(z,z_ijk(2),z_ijk(3));
r_ki = PointToEdgeDistance(z,z_ijk(3),z_ijk(1));

eps = 1e-20;
gamma_ij = zeros(size(r_ij));gamma_jk = gamma_ij;gamma_ki = gamma_ij;

gamma_ij(abs(z-z_ijk(1))<eps) = 1;
gamma_jk(abs(z-z_ijk(1))<eps) = 0;
gamma_ki(abs(z-z_ijk(1))<eps) = 0;

gamma_ij(abs(z-z_ijk(2))<eps) = 0;
gamma_jk(abs(z-z_ijk(2))<eps) = 1;
gamma_ki(abs(z-z_ijk(2))<eps) = 0;

gamma_ij(abs(z-z_ijk(3))<eps) = 0;
gamma_jk(abs(z-z_ijk(3))<eps) = 0;
gamma_ki(abs(z-z_ijk(3))<eps) = 1;

isnon0_denom_r = (r_jk.*r_ki+r_ij.*r_jk+r_ij.*r_ki)>eps;
denom_r_non0 = (r_jk(isnon0_denom_r).*r_ki(isnon0_denom_r)+...
    r_ij(isnon0_denom_r).*r_jk(isnon0_denom_r)+...
    r_ij(isnon0_denom_r).*r_ki(isnon0_denom_r));
gamma_ij(isnon0_denom_r) = (r_jk(isnon0_denom_r).*r_ki(isnon0_denom_r))./denom_r_non0;
gamma_jk(isnon0_denom_r) = (r_ij(isnon0_denom_r).*r_ki(isnon0_denom_r))./denom_r_non0;
gamma_ki(isnon0_denom_r) = (r_ij(isnon0_denom_r).*r_jk(isnon0_denom_r))./denom_r_non0;

% gamma_jk = (r_ij*r_ki)/(r_jk*r_ki+r_ij*r_jk+r_ij*r_ki);
% gamma_ki = (r_ij*r_jk)/(r_jk*r_ki+r_ij*r_jk+r_ij*r_ki);

% if abs(z-zi)<eps
%     gamma_ij = 1; gamma_jk = 0; gamma_ki = 0;
% elseif abs(z-zj)<eps
%     gamma_ij = 0; gamma_jk = 1; gamma_ki = 0;
% elseif abs(z-zk)<eps
%     gamma_ij = 0; gamma_jk = 0; gamma_ki = 1;
% else
%     gamma_ij = (r_jk*r_ki)/(r_jk*r_ki+r_ij*r_jk+r_ij*r_ki);
%     gamma_jk = (r_ij*r_ki)/(r_jk*r_ki+r_ij*r_jk+r_ij*r_ki);
%     gamma_ki = (r_ij*r_jk)/(r_jk*r_ki+r_ij*r_jk+r_ij*r_ki);
% end


end

function d = PointToEdgeDistance(z, z1, z2)
       n=size(z,1);
      pt = [real(z),imag(z),zeros(n,1)]; 
%       pt = [real(z),imag(z),0];       
      v1 = [real(z1),imag(z1),0];
      v2 = [real(z2),imag(z2),0];
%       a = v1 - v2;
        a = repmat(v1 - v2,n,1);
      b = pt - repmat(v2,n,1);
      d = normrow(cross(a,b)) ./ normrow(a);
end