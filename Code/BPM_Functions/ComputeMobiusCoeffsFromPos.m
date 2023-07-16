function [a,b,c,d] = ComputeMobiusCoeffsFromPos(z1,z2,z3,w1,w2,w3)
% function coeffs = ComputeMobiusCoeffsFromPos(z_pos,w_pos)
% z1 = z_pos(1);z2 = z_pos(2); z3 = z_pos(3);
% w1 = w_pos(1);w2 = w_pos(2); w3 = w_pos(3);
% Compute Mobius Coefficients From Positions
%   z1,z2,z3 are the original positins in the complex plane
%   w1,w2,w3 are the transformed positins in the complex plane

if or(isnan(z1),isnan(w1))
    a=nan; b=nan;c=nan; d=nan;
    return
end

a = det([z1*w1 w1 1; ...
         z2*w2 w2 1; ...
         z3*w3 w3 1]);
b = det([z1*w1 z1 w1; ...
         z2*w2 z2 w2; ...
         z3*w3 z3 w3]);
c = det([z1 w1 1; ...
         z2 w2 1; ...
         z3 w3 1]);
d = det([z1*w1 z1 1; ...
         z2*w2 z2 1; ...
         z3*w3 z3 1]);

norm_factor = sqrt(a*d-b*c);
if norm_factor < 1e-8
    norm_factor;
end
a = a/norm_factor; b = b/norm_factor; c = c/norm_factor; d = d/norm_factor;
% coeffs = [a,b,c,d];
end