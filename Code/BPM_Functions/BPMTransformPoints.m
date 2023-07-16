function [w,t_weights_ijkt,Mijk,Mij,Mjk,Mki] = BPMTransformPoints(Vz, Vw, F, Vz_subd, F_subd, orig_F_inds)
% function [w,t_weights_ijkt] = BPMTransformPoints(Vz, Vw, F, Vz_subd, F_subd, orig_F_inds)

%BPMTransformPoints Summary of this function goes here
%   Detailed explanation goes here

% extract the triangles star of each triangle
[F_star] = ExtractTrianglesStar(Vz, F, Vz_subd, F_subd);
[vf_ijk,vf_ij,vf_jk,vf_ki] = Extractvf(Vz,F,F_star);

w = zeros(size(Vz_subd,1),1);

t_weights_ijkt = zeros(size(Vz_subd,1),4);


Vzsubd_origF = zeros(size(Vz_subd,1),1); % for each v, one of the faces it belongs to
Vzsubd_origF(F_subd(:)) = orig_F_inds;

F_star1 = F_star(:,1); F_star2 = F_star(:,2); F_star3 = F_star(:,3);
wfij = nan(size(F_star,1),3);inds_fij = F(F_star1(~isnan(F_star1),1),:);
wfij(~isnan(F_star1),:) = [complex(Vw(inds_fij(:,1),1), Vw(inds_fij(:,1),2)),...
                           complex(Vw(inds_fij(:,2),1), Vw(inds_fij(:,2),2)),...
                           complex(Vw(inds_fij(:,3),1), Vw(inds_fij(:,3),2)),];
wfjk = nan(size(F_star,1),3);inds_fjk = F(F_star2(~isnan(F_star2),1),:);
wfjk(~isnan(F_star2),:) = [complex(Vw(inds_fjk(:,1),1), Vw(inds_fjk(:,1),2)),...
                           complex(Vw(inds_fjk(:,2),1), Vw(inds_fjk(:,2),2)),...
                           complex(Vw(inds_fjk(:,3),1), Vw(inds_fjk(:,3),2)),];
wfki = nan(size(F_star,1),3);inds_fki = F(F_star3(~isnan(F_star3),1),:);
wfki(~isnan(F_star3),:) = [complex(Vw(inds_fki(:,1),1), Vw(inds_fki(:,1),2)),...
                           complex(Vw(inds_fki(:,2),1), Vw(inds_fki(:,2),2)),...
                           complex(Vw(inds_fki(:,3),1), Vw(inds_fki(:,3),2)),];

for f=1:size(F_star,1)

inner_pts_ijk = Vz_subd(Vzsubd_origF==f,:);

[inner_pts_ijk_tr, vf_ijk_tr,vf_ij_tr,vf_jk_tr,vf_ki_tr] = ...
    IsometricEmbeddingFStar(inner_pts_ijk,vf_ijk(f,:),vf_ij(f,:),...
    vf_jk(f,:),vf_ki(f,:));

zf_ijk = [complex(vf_ijk_tr(:,1),vf_ijk_tr(:,2)),...
         complex(vf_ijk_tr(:,4),vf_ijk_tr(:,5)),...
         complex(vf_ijk_tr(:,7),vf_ijk_tr(:,8))];
zf_ij = [complex(vf_ij_tr(:,1),vf_ij_tr(:,2)),...
        complex(vf_ij_tr(:,4),vf_ij_tr(:,5)),...
        complex(vf_ij_tr(:,7),vf_ij_tr(:,8))];
zf_jk = [complex(vf_jk_tr(:,1),vf_jk_tr(:,2)),...
        complex(vf_jk_tr(:,4),vf_jk_tr(:,5)),...
        complex(vf_jk_tr(:,7),vf_jk_tr(:,8))];
zf_ki = [complex(vf_ki_tr(:,1),vf_ki_tr(:,2)),...
        complex(vf_ki_tr(:,4),vf_ki_tr(:,5)),...
        complex(vf_ki_tr(:,7),vf_ki_tr(:,8))];

wi = complex(Vw(F(f,1),1),Vw(F(f,1),2)); 
wj = complex(Vw(F(f,2),1),Vw(F(f,2),2));
wk = complex(Vw(F(f,3),1),Vw(F(f,3),2));


[a_ijk,b_ijk,c_ijk,d_ijk] = ComputeMobiusCoeffsFromPos(zf_ijk(1), zf_ijk(2), zf_ijk(3), wi, wj, wk);  



[a_ij,b_ij,c_ij,d_ij] = arrayfun(@ComputeMobiusCoeffsFromPos,  zf_ij(1), zf_ij(2), zf_ij(3), wfij(f,1), wfij(f,2), wfij(f,3));  

[a_jk,b_jk,c_jk,d_jk] = arrayfun(@ComputeMobiusCoeffsFromPos,  zf_jk(1), zf_jk(2), zf_jk(3), wfjk(f,1), wfjk(f,2), wfjk(f,3));  

[a_ki,b_ki,c_ki,d_ki] = arrayfun(@ComputeMobiusCoeffsFromPos,  zf_ki(1), zf_ki(2), zf_ki(3), wfki(f,1), wfki(f,2), wfki(f,3));  




Mijk = [a_ijk b_ijk; c_ijk d_ijk];
Mij = [a_ij b_ij; c_ij d_ij];
Mjk = [a_jk b_jk; c_jk d_jk];
Mki = [a_ki b_ki; c_ki d_ki];


z_inner_pts_tr = complex(inner_pts_ijk_tr(:,1),inner_pts_ijk_tr(:,2));
[O_mats,weights_ijkt] = BPMMatrixInterpolator(z_inner_pts_tr,zf_ijk,Mijk,Mij,Mjk,Mki);


w(Vzsubd_origF==f) = MobiusTransform(z_inner_pts_tr,O_mats(:,1),O_mats(:,2),O_mats(:,3),O_mats(:,4));
t_weights_ijkt(Vzsubd_origF==f,:) = weights_ijkt;
end

% t_weights_ijkt = 1-sum(t_weights_ijkt(:,1:3),2);
% t_weights_ijkt = t_weights_ijkt(:,3); 

% disp('******************************************')
% 
% % disply chosen weights
% for ind=[7,175,177]  
%     disp('*******************')
%     disp(['Vz coords:   ',num2str(Vz_subd(ind,[1,2]))])
%     disp(['   f=',num2str(Vzsubd_origF(ind)),'  coords:   ',num2str(vf_ijk(Vzsubd_origF(ind),[1,2,4,5,7,8]))])
%     disp([' ijkt weights:   ',num2str(t_weights_ijkt(ind,:))])
% end


end

function [vf_ijk,vf_ij,vf_jk,vf_ki] = Extractvf(Vz,F,F_star)
F_star1 = F_star(:,1); F_star2 = F_star(:,2); F_star3 = F_star(:,3);

% vertex positions of [ijk] 
inds_fijk = F;
vf_ijk = [Vz(inds_fijk(:,1),:) Vz(inds_fijk(:,2),:) Vz(inds_fijk(:,3),:)];

% vertex positions of the triangle adjacent [ij] with NaN where is no none
vf_ij = nan(size(F_star,1),9);inds_fij = F(F_star1(~isnan(F_star1),1),:);
vf_ij(~isnan(F_star1),:) = [Vz(inds_fij(:,1),:) Vz(inds_fij(:,2),:) Vz(inds_fij(:,3),:)];

% vertex positions of the triangle adjacent [jk] with NaN where is no none
vf_jk = nan(size(F_star,1),9);inds_fjk = F(F_star2(~isnan(F_star2),1),:);
vf_jk(~isnan(F_star2),:) = [Vz(inds_fjk(:,1),:) Vz(inds_fjk(:,2),:) Vz(inds_fjk(:,3),:)];

% vertex positions of the triangle adjacent [ki] with NaN where is no none
vf_ki = nan(size(F_star,1),9);inds_fki = F(F_star3(~isnan(F_star3),1),:);
vf_ki(~isnan(F_star3),:) = [Vz(inds_fki(:,1),:) Vz(inds_fki(:,2),:) Vz(inds_fki(:,3),:)];

end

function [inner_pts_ijk_tr, vf_ijk_tr,vf_ij_tr,vf_jk_tr,vf_ki_tr] = IsometricEmbeddingFStar(inner_pts_ijk,vf_ijk,vf_ij,vf_jk,vf_ki)
%  Inputs:
% inner_pts_ijk - nX3 interior points inside [ijk] to transform with [ijk] 
% vf_ijk - nX9 The vertex locations of the center triangle [ijk]
% vf_ij,f_jk,f_ki - 1X9 The vertex locations of the triangles adjacent to
%   edges ij,jk,ki
% vf_* is Nan if the center triangle is on the boundary and the
% corresponding triangle not exist
%  Outputs:
% vf_*_tr 1X9 the isometric embedding of the triangles star to the xy plane

% rotate the points to align [ijk] with xy plane
ijk_plane =  createPlane(vf_ijk(1:3), vf_ijk(4:6), vf_ijk(7:9));
TF = createBasisTransform3d('global', ijk_plane);
points = [vf_ijk(1:3);vf_ijk(4:6);vf_ijk(7:9);...
          vf_ij(1:3);vf_ij(4:6);vf_ij(7:9);...
          vf_jk(1:3);vf_jk(4:6);vf_jk(7:9);...
          vf_ki(1:3);vf_ki(4:6);vf_ki(7:9)]; 
points_tr1 = transformPoint3d(points, TF);    

inner_pts_ijk_tr = transformPoint3d(inner_pts_ijk, TF);


% rotate the neighbor triangles around edges ij, jk, ki resp to align them 
% with the XY plane

eij_dir = points_tr1(2,:) - points_tr1(1,:); 
ejk_dir = points_tr1(3,:) - points_tr1(2,:);
eki_dir = points_tr1(1,:) - points_tr1(3,:); 

if ~isnan(vf_ij)
    % transform [ijl] to align with ijk
    origin = points_tr1(1,:); %vi
    line = [origin eij_dir];
    n_ijl = cross(points_tr1(5,:)-points_tr1(4,:),points_tr1(6,:)-points_tr1(4,:));
    v1 = cross(eij_dir,[0 0 1]); v2 = n_ijl;
    theta = atan2(norm(cross(v1,v2)),dot(v1,v2))-pi/2;
    TF = createRotation3dLineAngle(line, theta);
    points_tr2 = points_tr1; 
    points_tr2(4:6,:) = transformPoint3d(points_tr1(4:6,:), TF);
else
    points_tr2 = points_tr1;
end
 if ~isnan(vf_jk)
    % transform [jkm] to align with ijk
    origin = points_tr2(2,:); %vj
    line = [origin ejk_dir];
    n_jkm = cross(points_tr2(8,:)-points_tr2(7,:),points_tr2(9,:)-points_tr2(7,:));
    v1 = cross(ejk_dir,[0 0 1]); v2 = n_jkm;
    theta = atan2(norm(cross(v1,v2)),dot(v1,v2))-pi/2;
    TF = createRotation3dLineAngle(line, theta);
    points_tr3 = points_tr2; 
    points_tr3(7:9,:) = transformPoint3d(points_tr2(7:9,:), TF);
else 
    points_tr3 = points_tr2;
end
if ~isnan(vf_ki)
    % transform [kin] to align with ijk
    origin = points_tr3(3,:); %vk
    line = [origin eki_dir];
    n_kin = cross(points_tr3(11,:)-points_tr3(10,:),points_tr3(12,:)-points_tr3(10,:));

    v1 = cross(eki_dir,[0 0 1]); v2 = n_kin;
    theta = atan2(norm(cross(v1,v2)),dot(v1,v2))-pi/2;
    TF = createRotation3dLineAngle(line, theta);
    points_tr4 = points_tr3; 
    points_tr4(10:12,:) = transformPoint3d(points_tr3(10:12,:), TF);
      
else 
    points_tr4 = points_tr3;
end

vf_ijk_tr = [points_tr4(1,:) points_tr4(2,:) points_tr4(3,:)];
vf_ij_tr = [points_tr4(4,:) points_tr4(5,:) points_tr4(6,:)];
vf_jk_tr = [points_tr4(7,:) points_tr4(8,:) points_tr4(9,:)];
vf_ki_tr = [points_tr4(10,:) points_tr4(11,:) points_tr4(12,:)];






end

function [F_star] = ExtractTrianglesStar(Vz, F, Vz_subd, F_subd)
% extract the triangles star of each triangle
% Inputs:
% Outputs:
% F_star- [f_ijk, f_ij, f_jk, f_ki] by this order
% F_star(i,j) = NaN if there is no neighbor to the corresponding edge
% Vsubd_origF1 -  
 % compute neighbors of each face ijk to extract Mijk,Mij,Mjk,Mki 
TR = triangulation(F,Vz);
% compute F_star
F_star = neighbors(TR); % ordered by ijk,jk,ki,ij
F_star = circshift(F_star,1,2); % order by ij, jk, ki

end



% Todo:
% take the vertices from faces F_subd to track to which face each vertex
% belongs at the cost of computing 3 times for each vertex (can be 
% optimized later using unique...)