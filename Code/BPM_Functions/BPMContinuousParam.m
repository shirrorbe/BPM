function [Mw_fine, Mz_fine,t_weights_ijkt,Mijk,Mij,Mjk,Mki] = BPMContinuousParam(Vz, Vw, F ,n_subd)
% function [Mw_fine, Mz_fine,t_weights_ijkt] = BPMContinuousParam(Vz, Vw, F ,n_subd)

% sampling the domain by subdividing the mesh
if n_subd > 0
    [Vz_subd, F_subd] = meshSubdivision(Vz, F);
    for i=1:(n_subd-1)
        [Vz_subd, F_subd] = meshSubdivision(Vz_subd, F_subd);
    end
    % orig_F : the original face that each new face lies in
    orig_F = repmat((1:size(F,1)),4^n_subd,1);
    orig_F = orig_F(:);
    orig_F = repmat(orig_F,3,1); % the original face that each new face lies in
else
    Vz_subd = Vz; F_subd=F;
    orig_F = (1:size(F,1))';
end
Mz_fine.V = Vz_subd; Mz_fine.F = F_subd;


 % compute neighbors of each face ijk to extract Mijk,Mij,Mjk,Mki 
% TR = triangulation(F,Vz);
% % compute F_star
% % F_star- [f_ijk, f_ij, f_jk, f_ki] by this order
% % F_star(i,j) = NaN if there is no neighbor to the corresponding edge
% F_star = neighbors(TR); % ordered by ijk,jk,ki,ij
% F_star = circshift(F_star,1,2); % order by ij, jk, ki
% % take the vertices from faces F_subd to track to which face each vertex
% % belongs at the cost of computing 3 times for each vertex (can be 
% % optimized later using unique...)
% Vsubd_origF = Vz_subd([F_subd(:,1);F_subd(:,2);F_subd(:,3)],:);

% [w,t_weights_ijkt] = BPMTransformPoints(Vz, Vw, F, Vz_subd, F_subd, orig_F);
[w,t_weights_ijkt,Mijk,Mij,Mjk,Mki] = BPMTransformPoints(Vz, Vw, F, Vz_subd, F_subd, orig_F);

% [w] = GeneralizedPCMCoeffs(F_star, F, Vsubd_origF, orig_F,Vz, Vw, F_subd);
% Vw_fine=zeros(size(Vz_subd));
% Vw_fine(F_subd(:),:) = [real(w) imag(w) zeros(size(w,1),1)];
Vw_fine = [real(w) imag(w) zeros(size(w,1),1)];
Mw_fine.V = Vw_fine; Mw_fine.F = F_subd;

end