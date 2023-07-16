function [w] = LinearInterp(z ,zi, zj, zk, wi, wj, wk)
% LinearInterp
% Inputs:
%  z- the points to tranform
%  zi, zj, zk - the end points of the middle triangle.
% Output:
% w - the linear interpolation inside triangle ijk 

v = [real(z) imag(z)];
vi = repmat([real(zi) imag(zi)],size(z,1),1);
vj = repmat([real(zj) imag(zj)],size(z,1),1);
vk = repmat([real(zk) imag(zk)],size(z,1),1);
B = barycentric_coordinates(v,vi,vj,vk);
w = (B(:,1)*wi + B(:,2)*wj + B(:,3)*wk)./(B(:,1) + B(:,2) + B(:,3));
end
