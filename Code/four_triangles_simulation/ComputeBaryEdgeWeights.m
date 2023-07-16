function [gamma_ij,gamma_jk,gamma_ki] = ComputeBaryEdgeWeights(Mz,Mw,vi,vj,vk,n_subd)
Vz = Mz.V; Vw = Mw.V; F=Mz.F;
[Mw_fine, Mz_fine] = BPMContinuousParam(Vz, Vw, F ,n_subd);
Vz_fine = Mz_fine.V; Vw_fine = Mw_fine.V; F_fine = Mz_fine.F;

% a =  repmat(vi - vj,size(Vz_fine,1),1); b = Vz_fine - vj;
% d = normrow(cross(a,b)) ./ normrow(a);
% d = d/max(d);

% w_ij = zeros(size(d));
% w_ij(d>0) = 1./d(d>0);
% w_ij = w_ij/max(w_ij);
% w_ij(d<=0) = 1;
% 
% w_jk = zeros(size(d));
% w_jk(d>0) = 1./d(d>0);
% w_jk = w_jk/max(w_jk);
% w_jk(d<=0) = 1;
% 
% w_ki = zeros(size(d));
% w_ki(d>0) = 1./d(d>0);
% w_ki = w_ki/max(w_ki);
% w_ki(d<=0) = 1;
% 
% 
% figure
% title('weights eij');
% nf_fine = size(F_fine,1); nff1 = floor(nf_fine/4);
% hp = patch('Faces',F_fine(1:nff1,:),'Vertices',Vz_fine(:,1:2),'FaceColor','w'); axis equal; hold on; 
% set(hp,'FaceVertexCData',w_ij,'FaceColor','interp');
% % draw_point_2d([Vz(:,1),Vz(:,2)],'MarkerSize',30);
% % draw_point_2d([Vz_fine(d<=0,1),Vz_fine(d<=0,2)],'MarkerSize',30);
% 
% aa = axis; aa = aa + 0.5*[-1,1,-1,1]; axis(aa);
% set(gcf,'WindowStyle','docked')
% colorbar
% 

eps = 1e-20;
r_ij = PointToEdgeDistance(Vz_fine, vi, vj);
r_jk = PointToEdgeDistance(Vz_fine, vj, vk);
r_ki = PointToEdgeDistance(Vz_fine, vk, vi);

% gamma_ij = zeros(size(r_ij));
% gamma_ij(normrow(Vz_fine-vi)>eps) = 1./r_ij(r_ij>eps);

% gamma_ij(r_ij>eps) = 1./r_ij(r_ij>eps);

% gamma_jk = zeros(size(r_jk));
% gamma_jk(r_jk>eps) = 1./r_jk(r_jk>eps);

% gamma_ki = zeros(size(r_ki));
% gamma_ki(r_ki>eps) = 1./r_ki(r_ki>eps);

gamma_ij = (r_jk.*r_ki)./(r_jk.*r_ki+r_ij.*r_jk+r_ij.*r_ki);
gamma_jk = (r_ij.*r_ki)./(r_jk.*r_ki+r_ij.*r_jk+r_ij.*r_ki);
gamma_ki = (r_ij.*r_jk)./(r_jk.*r_ki+r_ij.*r_jk+r_ij.*r_ki);


gamma_ij(normrow(Vz_fine-vi)<=eps) = 1;gamma_jk(normrow(Vz_fine-vi)<=eps) = 0;gamma_ki(normrow(Vz_fine-vi)<=eps) = 0;
gamma_ij(normrow(Vz_fine-vj)<=eps) = 0;gamma_jk(normrow(Vz_fine-vj)<=eps) = 1;gamma_ki(normrow(Vz_fine-vj)<=eps) = 0;
gamma_ij(normrow(Vz_fine-vk)<=eps) = 0;gamma_jk(normrow(Vz_fine-vk)<=eps) = 0;gamma_ki(normrow(Vz_fine-vk)<=eps) = 1;




% r_ij = PointToEdgeDistance(Vz_fine, vi, vj);
% gamma_ij = zeros(size(r_ij));
% % gamma_ij(Vz_fine->eps) = 1./r_ij(r_ij>eps);
% gamma_ij(r_ij>eps) = 1./r_ij(r_ij>eps);
% 
% r_jk = PointToEdgeDistance(Vz_fine, vj, vk);
% gamma_jk = zeros(size(r_jk));
% gamma_jk(r_jk>eps) = 1./r_jk(r_jk>eps);
% 
% r_ki = PointToEdgeDistance(Vz_fine, vk, vi);
% gamma_ki = zeros(size(r_ki));
% gamma_ki(r_ki>eps) = 1./r_ki(r_ki>eps);
% 
% gamma_ij(r_ij<=eps) = 1;gamma_jk(r_ij<=eps) = 0;gamma_ki(r_ij<=eps) = 0;
% gamma_ij(r_jk<=eps) = 0;gamma_jk(r_jk<=eps) = 1;gamma_ki(r_jk<=eps) = 0;
% gamma_ij(r_ki<=eps) = 0;gamma_jk(r_ki<=eps) = 0;gamma_ki(r_ki<=eps) = 1;
% 
% sum_gamma = gamma_ij+gamma_jk+gamma_ki;
% gamma_ij = gamma_ij./sum_gamma;
% gamma_jk = gamma_jk./sum_gamma;
% gamma_ki = gamma_ki./sum_gamma;



% 
% 
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
% 
% 






end



function d = PointToEdgeDistance(pnt, vi, vj)
Vz_fine = pnt;
a =  repmat(vi - vj,size(Vz_fine,1),1); b = Vz_fine - vj;
d = normrow(cross(a,b)) ./ normrow(a);
end