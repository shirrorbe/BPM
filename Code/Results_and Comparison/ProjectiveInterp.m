 function [w] = ProjectiveInterp(z ,zi, zj, zk, wi, wj, wk, ui, uj, uk)
% ProjectiveInterp
% Inputs:
%  z- the points to tranform
%  zi, zj, zk - the end points of the middle triangle.
% Output:
% w - the blended Mobius transformation of the z points

v = [real(z) imag(z)];
vi = repmat([real(zi) imag(zi)],size(z,1),1);
vj = repmat([real(zj) imag(zj)],size(z,1),1);
vk = repmat([real(zk) imag(zk)],size(z,1),1);
B = barycentric_coordinates(v,vi,vj,vk);
w = (B(:,1)*exp(-ui)*wi + B(:,2)*exp(-uj)*wj + B(:,3)*exp(-uk)*wk)./(B(:,1)*exp(-ui) + B(:,2)*exp(-uj) + B(:,3)*exp(-uk));
end
