function  h=function_visualization_on_mesh(v,f,f_vals)
if size(v,1)==size(f_vals,1)
    f_clr='interp';
else
    f_clr='flat';
end

h=patch('Faces',f,'Vertices',v,'FaceVertexCData',full(f_vals),'FaceColor',f_clr,'EdgeColor','none');
axis vis3d
colorbar
axis equal
end

