function [f_areas,v_areas] = calculate_face_vertex_areas  (v,f)

nv=size(v,1);
nf=size(f,1);


cross_prod = cross(v(f(:,2),:)-v(f(:,1),:),v(f(:,3),:)-v(f(:,1),:),2);
f_areas = 0.5*sqrt(sum(cross_prod.^2,2));
vf_adj = sparse(f, repmat((1:nf)', 1, 3), repmat(f_areas, 1, 3), nv, nf);

v_areas=sum(vf_adj,2)/3;

end

