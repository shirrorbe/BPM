% generate an obj with texture coordinates M_vt of MeshWrap object M.
function generate_texture_mesh(M, M_vt, texture_im_name,...
    out_obj_name)

options.object_texture = M_vt;
options.nm_file = texture_im_name;

path = '';
if any(out_obj_name=='\')
    path = out_obj_name(1:find(out_obj_name=='\', 1, 'last'));
    out_obj_name = out_obj_name(find(out_obj_name=='\', 1, 'last') + 1:end);
end
% options.vertex_normals = M.vertex_normals();
options.vertex_normals = [];

write_obj2(path, out_obj_name, M.vertices, M.triangles, options);
