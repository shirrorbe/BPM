function [gamma_ij,gamma_jk,gamma_ki,r_ij,r_jk,r_ki] = BarycentricWeights(z,zi,zj,zk)


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
end


function d = PointToEdgeDistance(z, z1, z2)
      pt = [real(z),imag(z),0];       
      v1 = [real(z1),imag(z1),0];
      v2 = [real(z2),imag(z2),0];
      a = v1 - v2;
      b = pt - v2;
      d = norm(cross(a,b)) / norm(a);
end

