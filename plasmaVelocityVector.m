function [v,cov] = plasmaVelocityVector(vlos,vloserr,k)
%
% [v,cov] = plasmaVelocityVector(vlos,vloserr,k)
%
% Calculate plasma drift velocity vector from line-of-sight
% velocity estimates.
%
% INPUT:
%  vlos    a column vector of line-of-sight velocity estimates (m/s)
%  vloserr standard deviations of vlos (m/s)
%  k       scattering wave vector directions (unit vectors) (n x 3
%          matrix) in a right-handed cartesian coordinate system
%
% OUTPUT:
%  v       1 x 3 velocity vector (m/s)
%  cov     3 x 3 error covariance matrix of v (m^s^-2)
%
%
% IV 2016
%

% Since each row of k is a unit vector in direction of scattering
% wave vector, we have vlos = kv + e, where std(e) = vloserr.
% Thus

prinf = 1./(vloserr.^2);
try
    % try to invert the posterior information matrix
    cov = inv(k' * ( diag(prinf) * k));
catch exception
    % if the inversion did not work, we will just give up...
    cov = NaN(3);
end
v = cov * k' * (vlos.*prinf);

end