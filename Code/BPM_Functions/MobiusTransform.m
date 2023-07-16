function [Z_m, U] = MobiusTransform(Z,a,b,c,d)
%Compute the Mobius transform for every z in X (z is a complex number)
% a,b,c,d are normalized such that a*d-b*c=1
% Z, Z_m is in C^(mXn)
% Z_m = (az+b)(cz+d)^-1

% output:
% U the scaling factors

Z_m = (a.*Z+b)./(c.*Z+d);
U = (abs(c.*Z+d)).^-4;
end