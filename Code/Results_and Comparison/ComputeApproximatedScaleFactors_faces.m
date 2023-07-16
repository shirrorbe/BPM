function u = ComputeApproximatedScaleFactors_faces(Vs,Vt,F)
% approximate the scale factors by using the notion to discrete conformal equivalence
% of the paper "Conformal Equivalence Of Triangle Meshes" by Springborn et
% al.
% for each triangle e^u_i = (l_ij_tilde/l_ij) * (l_kj/l_kj_tilde) * (l_ik_tilde/l_ik)

l_ij_tilde = normrow(Vs(F(:,1),:)-Vs(F(:,2),:));
l_ik_tilde = normrow(Vs(F(:,1),:)-Vs(F(:,3),:));
l_kj_tilde = normrow(Vs(F(:,3),:)-Vs(F(:,2),:));

l_ij = normrow(Vt(F(:,1),:)-Vt(F(:,2),:));
l_ik = normrow(Vt(F(:,1),:)-Vt(F(:,3),:));
l_kj = normrow(Vt(F(:,3),:)-Vt(F(:,2),:));

% l_ij_tilde = normrow(Vt(F(:,1),:)-Vt(F(:,2),:));
% l_ik_tilde = normrow(Vt(F(:,1),:)-Vt(F(:,3),:));
% l_kj_tilde = normrow(Vt(F(:,3),:)-Vt(F(:,2),:));
% 
% l_ij = normrow(Vs(F(:,1),:)-Vs(F(:,2),:));
% l_ik = normrow(Vs(F(:,1),:)-Vs(F(:,3),:));
% l_kj = normrow(Vs(F(:,3),:)-Vs(F(:,2),:));
% 

log_ij = log(l_ij_tilde./l_ij);
log_ik = log(l_ik_tilde./l_ik);
log_kj = log(l_kj_tilde./l_kj);

ui = log_ij - log_kj + log_ik;
uj = log_ij + log_kj - log_ik;
uk = - log_ij + log_kj + log_ik;

u = [ui uj uk];


end

