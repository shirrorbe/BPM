function [w] = GeneralizedPCMCoeffs(F_star, F, Vz1_faces, orig_F,Vz, Vw,F1)
% F_star - the neighbors of each triangle ordered by the edges jk, ki, ij
% Vz1_faces - vertices of each face in F1 contaminated by the 3 columns
% compute z for each face star
% flatten each face star to the complex plane

F_star1 = F_star(:,1); F_star2 = F_star(:,2); F_star3 = F_star(:,3);

% vertex positions of [ijk] 
inds_fijk = F;
vfijk = [Vz(inds_fijk(:,1),:) Vz(inds_fijk(:,2),:) Vz(inds_fijk(:,3),:)];

% vertex positions of the triangle adjacent [ij] with NaN where is no none
vfij = nan(size(F_star,1),9);inds_fij = F(F_star1(~isnan(F_star1),1),:);
vfij(~isnan(F_star1),:) = [Vz(inds_fij(:,1),:) Vz(inds_fij(:,2),:) Vz(inds_fij(:,3),:)];

% vertex positions of the triangle adjacent [jk] with NaN where is no none
vfjk = nan(size(F_star,1),9);inds_fjk = F(F_star2(~isnan(F_star2),1),:);
vfjk(~isnan(F_star2),:) = [Vz(inds_fjk(:,1),:) Vz(inds_fjk(:,2),:) Vz(inds_fjk(:,3),:)];

% vertex positions of the triangle adjacent [ki] with NaN where is no none
vfki = nan(size(F_star,1),9);inds_fki = F(F_star3(~isnan(F_star3),1),:);
vfki(~isnan(F_star3),:) = [Vz(inds_fki(:,1),:) Vz(inds_fki(:,2),:) Vz(inds_fki(:,3),:)];

% Isometric enbedding of [ijk] with its neighbors on edges [ij],[jk],[ki]
[zfijk,zfij, zfjk, zfki, z1] = UnfoldTrStarToComplexPlane(vfijk,vfij, vfjk, vfki, Vz1_faces, orig_F);


wi = complex(Vw(F(:,1),1),Vw(F(:,1),2)); 
wj = complex(Vw(F(:,2),1),Vw(F(:,2),2));
wk = complex(Vw(F(:,3),1),Vw(F(:,3),2));


[a_ijk,b_ijk,c_ijk,d_ijk] = arrayfun(@ComputeMobiusCoeffsFromPos, zfijk(:,1), zfijk(:,2), zfijk(:,3), wi, wj, wk);  



wfij = nan(size(F_star,1),3);inds_fij = F(F_star1(~isnan(F_star1),1),:);
wfij(~isnan(F_star1),:) = [complex(Vw(inds_fij(:,1),1), Vw(inds_fij(:,1),2)),...
                           complex(Vw(inds_fij(:,2),1), Vw(inds_fij(:,2),2)),...
                           complex(Vw(inds_fij(:,3),1), Vw(inds_fij(:,3),2)),];
[a_ij,b_ij,c_ij,d_ij] = arrayfun(@ComputeMobiusCoeffsFromPos,  zfij(:,1), zfij(:,2), zfij(:,3), wfij(:,1), wfij(:,2), wfij(:,3));  

wfjk = nan(size(F_star,1),3);inds_fjk = F(F_star2(~isnan(F_star2),1),:);
wfjk(~isnan(F_star2),:) = [complex(Vw(inds_fjk(:,1),1), Vw(inds_fjk(:,1),2)),...
                           complex(Vw(inds_fjk(:,2),1), Vw(inds_fjk(:,2),2)),...
                           complex(Vw(inds_fjk(:,3),1), Vw(inds_fjk(:,3),2)),];
[a_jk,b_jk,c_jk,d_jk] = arrayfun(@ComputeMobiusCoeffsFromPos,  zfjk(:,1), zfjk(:,2), zfjk(:,3), wfjk(:,1), wfjk(:,2), wfjk(:,3));  


wfki = nan(size(F_star,1),3);inds_fki = F(F_star3(~isnan(F_star3),1),:);
wfki(~isnan(F_star3),:) = [complex(Vw(inds_fki(:,1),1), Vw(inds_fki(:,1),2)),...
                           complex(Vw(inds_fki(:,2),1), Vw(inds_fki(:,2),2)),...
                           complex(Vw(inds_fki(:,3),1), Vw(inds_fki(:,3),2)),];
[a_ki,b_ki,c_ki,d_ki] = arrayfun(@ComputeMobiusCoeffsFromPos,  zfki(:,1), zfki(:,2), zfki(:,3), wfki(:,1), wfki(:,2), wfki(:,3));  



a1_ijk = a_ijk(orig_F); b1_ijk = b_ijk(orig_F);c1_ijk = c_ijk(orig_F); d1_ijk = d_ijk(orig_F);

a1_ij = a_ij(orig_F); b1_ij = b_ij(orig_F);c1_ij = c_ij(orig_F); d1_ij = d_ij(orig_F);

a1_jk = a_jk(orig_F); b1_jk = b_jk(orig_F);c1_jk = c_jk(orig_F); d1_jk = d_jk(orig_F);

a1_ki = a_ki(orig_F); b1_ki = b_ki(orig_F);c1_ki = c_ki(orig_F); d1_ki = d_ki(orig_F);
   
zi= zfijk(:,1); zj= zfijk(:,2); zk= zfijk(:,3);
z1_i = zi(orig_F); z1_j = zj(orig_F); z1_k = zk(orig_F);

% [w] = arrayfun(@SlerpBlendMobius, z1 ,z1_i, z1_j, z1_k, a1_ijk,b1_ijk,c1_ijk,d1_ijk,...
%         a1_ij,b1_ij,c1_ij,d1_ij, a1_jk,b1_jk,c1_jk,d1_jk, a1_ki,b1_ki,c1_ki,d1_ki);

[w] = arrayfun(@SlerpBlendMobius_stiff, z1 ,z1_i, z1_j, z1_k, a1_ijk,b1_ijk,c1_ijk,d1_ijk,...
        a1_ij,b1_ij,c1_ij,d1_ij, a1_jk,b1_jk,c1_jk,d1_jk, a1_ki,b1_ki,c1_ki,d1_ki);



end




function [zfijk, zfij, zfjk, zfki, z1] = UnfoldTrStarToComplexPlane(vfijk,vfij, vfjk, vfki, Vz1_faces,orig_F)
% vf* - n by 9 triangles vertices [v1 v2 v3] 
n_pts = size(vfijk,1);

% treat nan ([ijk] is on boundary and one of its triangle stars don't exist)
vfij_old = vfij; vfjk_old = vfjk; vfki_old = vfki;
% vfij(isnan(vfij)) = 0;vfjk(isnan(vfjk)) = 0;vfki(isnan(vfki)) = 0; % dummy 


% rotate the points of the 4 triangles such that [ijk] is allgined with the XY plane

vfijk_tr = zeros(size(vfijk)); vfij_tr = zeros(size(vfij));
vfjk_tr = zeros(size(vfjk)); vfki_tr = zeros(size(vfki));

Vz1_faces_tr = zeros(size(Vz1_faces));
for i=1:n_pts
    % rotate the points to align [ijk] with xy plane
    ijk_plane =  createPlane(vfijk(i,1:3), vfijk(i,4:6), vfijk(i,7:9));
    TF = createBasisTransform3d('global', ijk_plane);
    points = [vfijk(i,1:3);vfijk(i,4:6);vfijk(i,7:9);...
              vfij(i,1:3);vfij(i,4:6);vfij(i,7:9);...
              vfjk(i,1:3);vfjk(i,4:6);vfjk(i,7:9);...
              vfki(i,1:3);vfki(i,4:6);vfki(i,7:9)]; 
    newPoints = transformPoint3d(points, TF);    
    new_interior_pts = transformPoint3d(Vz1_faces(orig_F==i,:), TF);
    Vz1_faces_tr(orig_F==i,:) = new_interior_pts;
    

    % rotate the neighbor triangles around edges ij, jk, ki resp. to flatten
    % the triangles in each star to the xy plane
    eij_dir = newPoints(2,:) - newPoints(1,:); %eij
    ejk_dir = newPoints(3,:) - newPoints(2,:); %ejk
    eki_dir = newPoints(1,:) - newPoints(3,:); %ejk

    if ~isnan(vfij(i,:))
        % transform [ijl] to align with ijk
        origin = newPoints(1,:); %vi
%         direction = newPoints(2,:) - newPoints(1,:); %eij
%         direction = -direction;
        line = [origin eij_dir];
       
        n_ijl = cross(newPoints(5,:)-newPoints(4,:),newPoints(6,:)-newPoints(4,:));

        v1 = cross(eij_dir,[0 0 1]); v2 = n_ijl;
        theta = atan2(norm(cross(v1,v2)),dot(v1,v2))-pi/2;

% %         acos(dot(v1 / norm(v1), v2 / norm(v2)))
%         if theta<0
%         end

        TF = createRotation3dLineAngle(line, theta);
        newPoints2 = newPoints; 
        newPoints2(4:6,:) = transformPoint3d(newPoints(4:6,:), TF);
    else
        newPoints2 = newPoints;
    end
    
    if ~isnan(vfjk(i,:))
        % transform [jkm] to align with ijk
        origin = newPoints2(2,:); %vj
%         direction = newPoints2(3,:) - newPoints2(2,:); %ejk
%         direction = -direction;
%         line = [origin direction];
        line = [origin ejk_dir];
        % jk_plane = createPlane(newPoints2(7,:), newPoints2(8,:), newPoints2(9,:));
        % n_jkm = cross(jk_plane(4:6),jk_plane(7:9));
        n_jkm = cross(newPoints2(8,:)-newPoints2(7,:),newPoints2(9,:)-newPoints2(7,:));
%         % n_jkm = cross(newPoints2(9,:)-newPoints2(9,:)-,jk_plane(7:9))
%         theta = atan2(norm(cross(n_ijk,n_jkm)),dot(n_ijk,n_jkm));
% %         theta = -acos(dot(n_ijk,n_jkm)./(norm(n_ijk)*norm(n_jkm)));
        v1 = cross(ejk_dir,[0 0 1]); v2 = n_jkm;
        theta = atan2(norm(cross(v1,v2)),dot(v1,v2))-pi/2;
        TF = createRotation3dLineAngle(line, theta);
        newPoints3 = newPoints2; 
        newPoints3(7:9,:) = transformPoint3d(newPoints2(7:9,:), TF);
    else 
        newPoints3 = newPoints2;
    end

    if ~isnan(vfki(i,:))
        % transform [kin] to align with ijk
        origin = newPoints3(3,:); %vk
%         direction = newPoints3(1,:) - newPoints3(3,:); %eik
%         direction = -direction;
%         line = [origin direction];
        line = [origin eki_dir];

%         n_ijk = cross(newPoints3(2,:) - newPoints3(1,:),newPoints3(3,:) - newPoints3(1,:));
        % ki_plane = createPlane(newPoints3(10,:), newPoints3(11,:), newPoints3(12,:));
        % n_kin = cross(ki_plane(4:6),ki_plane(7:9));
        n_kin = cross(newPoints3(11,:)-newPoints3(10,:),newPoints3(12,:)-newPoints3(10,:));
%         theta = atan2(norm(cross(n_ijk,n_kin)),dot(n_ijk,n_kin));
%         theta = acos(dot(n_ijk,n_kin)./(norm(n_ijk)*norm(n_kin)));

        v1 = cross(eki_dir,[0 0 1]); v2 = n_kin;
        theta = atan2(norm(cross(v1,v2)),dot(v1,v2))-pi/2;
        TF = createRotation3dLineAngle(line, theta);
        newPoints4 = newPoints3; 
        newPoints4(10:12,:) = transformPoint3d(newPoints3(10:12,:), TF);
          
    else 
        newPoints4 = newPoints3;
    end
 
    pts_tr = newPoints4;
    vfijk_tr(i,:) = [pts_tr(1,:) pts_tr(2,:) pts_tr(3,:)];
    vfij_tr(i,:) = [pts_tr(4,:) pts_tr(5,:) pts_tr(6,:)];
    vfjk_tr(i,:) = [pts_tr(7,:) pts_tr(8,:) pts_tr(9,:)];
    vfki_tr(i,:) = [pts_tr(10,:) pts_tr(11,:) pts_tr(12,:)];
end



zfijk = [complex(vfijk_tr(:,1),vfijk_tr(:,2)),...
         complex(vfijk_tr(:,4),vfijk_tr(:,5)),...
         complex(vfijk_tr(:,7),vfijk_tr(:,8))];
zfij = [complex(vfij_tr(:,1),vfij_tr(:,2)),...
        complex(vfij_tr(:,4),vfij_tr(:,5)),...
        complex(vfij_tr(:,7),vfij_tr(:,8))];
zfjk = [complex(vfjk_tr(:,1),vfjk_tr(:,2)),...
        complex(vfjk_tr(:,4),vfjk_tr(:,5)),...
        complex(vfjk_tr(:,7),vfjk_tr(:,8))];
zfki = [complex(vfki_tr(:,1),vfki_tr(:,2)),...
        complex(vfki_tr(:,4),vfki_tr(:,5)),...
        complex(vfki_tr(:,7),vfki_tr(:,8))];

% treat nan ([ijk] is on boundary and one of its triangle stars don't exist)
% delete dummy values
% zfij(isnan(vfij_old)) = nan;
isnan_ij = isnan(vfij_old(:,1));
zfij(isnan_ij,:) = nan(sum(isnan_ij),3);

isnan_jk = isnan(vfjk_old(:,1));
zfjk(isnan_jk,:) = nan(sum(isnan_jk),3);

isnan_ki = isnan(vfki_old(:,1));
zfki(isnan_ki,:) = nan(sum(isnan_ki),3);

 z1 = complex(Vz1_faces_tr(:,1),Vz1_faces_tr(:,2));


end




