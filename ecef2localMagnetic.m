function [kBned] = ecef2localMagnetic( k , ecefS , llhgS , Bned )
%
% [kB] = ecef2localMagnetic( k , ecefS , llhgS , Bned )
%
% convert k-vectors from ecef coordinates to *magnetic* enu
%
% INPUT:
%  k      the k-vectors in cartesian ecef coordinates, as row vectors
%  ecefS  ecef coordinates of the scattering volume
%  llhgS  geodetic latitude, longitude, height of the scattering
%         volume
%  Bned   magnetic field vector (nT) in local ned coordinates
%
% OUTPUT:
%  kBned  the k-vectors in local magnetic north-east-down coordinates
%
%
% IV 2017
%

    if nargin < 5
        show = false;
    end

% allocate variables...
kned = zeros(size(k));
kBned = zeros(size(k));

% number of k-vectors
nk = size(k,1);

% the wgs84 ellipsoid
wgs84 = wgs84Ellipsoid('meters');

% geodetic ned coordinates
for ik = 1:nk

    % this conversion is numerically unstable since we add about one meter to thousands of kilometers... should multiply k with a large number!
    X = ecefS(ik,1) + k(ik,1)*1e4;
    Y = ecefS(ik,2) + k(ik,2)*1e4;
    Z = ecefS(ik,3) + k(ik,3)*1e4;
    [n,e,d] = ecef2ned(X, Y, Z, llhgS(ik,1), llhgS(ik,2), llhgS(ik,3), wgs84 ...
                       );
    kned(ik,:) = [n,e,d]/1e4;

    % make sure that this is a unit vector...
    kned(ik,:) = kned(ik,:)./sqrt(sum(kned(ik,:).^2));
end

% geomagnetic ned (magneticn north, magnetic east, along B)
for ik = 1:nk
    % horizontal component of magnetic north
    Bnh = [Bned(ik,1) Bned(ik,2) 0 ];
    % magnetic down....
    Bd = Bned(ik,:);
    % magnetic east (Bd x Bnh)
    Be = cross( Bd , Bnh );
    % magnetic north
    Bn = cross( Be , Bd );
    % normalize all to unit length
    Bn = Bn./sqrt(sum(Bn.^2));
    Be = Be./sqrt(sum(Be.^2));
    Bd = Bd./sqrt(sum(Bd.^2));

    % project the k-vector to the magnetic coordinates
    kBned(ik,1) = sum(kned(ik,:).*Bn);
    kBned(ik,2) = sum(kned(ik,:).*Be);
    kBned(ik,3) = sum(kned(ik,:).*Bd);
end

end
