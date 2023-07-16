function [Vz,Vw,F] = ExtractMeshesFromUV(filename)
obj = readObj([filename '.obj']);
UV = zeros(size(obj.v,1),2);
UV(obj.f.v,:) = obj.vt(obj.f.vt,1:2);
% UV(obj.f.v,:) = obj.vt(obj.f.vt,:);

Mz0.V = obj.v; Mz0.F = obj.f.v;
Mw0.V = [UV(:,1:2), zeros(size(UV,1),1)]; Mw0.F = Mz0.F;
% out_name=filename;
% Tri = Mz0.F'; Pts = Mz0.V'; out = [out_name,'_out.off'];
% data_to_off(Tri,Pts,out)
% Tri = Mw0.F'; Pts = Mw0.V'; out = [out_name,'_in.off'];
% data_to_off(Tri,Pts,out)
Vz = Mz0.V; Vw = Mw0.V; F = Mz0.F;
end

