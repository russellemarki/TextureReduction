function [ori_out] = ori2eulerFZ(ori)
% finds all crystallographically equivalent orientations and then selects
% those in the Euler fundamental zone
%
% Syntax
%   [phi] = ori2eulerFZ(ori)    % Bunge fundamental region
%   [abg] = ori2eulerFZ(ori, 'Matthies')  % Matthies fundamental region
%
% Input
%  ori        - @orientation
%  CS1, CS2 - crystal @symmetry
%
% Output
%  phi - [phi1, Phi, phi2] - Euler angles
%

CS = ori.CS;
SS = ori.SS;

ori = ori(:);

% get the fundamental region
[maxphi1,maxPhi,maxphi2] = fundamentalRegionEuler(CS, SS);

% symmetrise
ori_sym = (SS.rotation_special * quaternion(ori)).' * CS.rotation_special;

ori_sym = reshape(ori_sym, length(ori), []);

[phi1_sym, Phi_sym, phi2_sym] = Euler(ori_sym);

% the simple part
phi1_sym = mod(phi1_sym,2*pi/SS.multiplicityZ);
phi2_sym = mod(phi2_sym,2*pi/CS.multiplicityZ);

mult = 1 - 1e-12;

if CS.id == 45
    %special check for cubic symmetry, roof fundamental zone
    isInside = (phi1_sym <= maxphi1) & (...
        (phi2_sym <= pi/4 & Phi_sym >= acos(sin(phi2_sym)./sqrt(1 + sin(phi2_sym).^2)) & Phi_sym <= pi/2) | ...
        (phi2_sym >  pi/4 & Phi_sym >= acos(cos(phi2_sym)./sqrt(1 + cos(phi2_sym).^2)) & Phi_sym <= pi/2 & phi2_sym <= pi/2));
    
    
    %if there aren't any matches we relax the constraints a little bit
    nind = find(~any(isInside, 2));
    
    isInside(nind, :) = (phi1_sym(nind, :) <= maxphi1) & (...
        (phi2_sym(nind, :) <= pi/4 & Phi_sym(nind, :) >= acos(sin(phi2_sym(nind, :))./sqrt(1 + sin(phi2_sym(nind, :)).^2))*mult & Phi_sym(nind, :) <= pi/2) | ...
        (phi2_sym(nind, :) >  pi/4 & Phi_sym(nind, :) >= acos(cos(phi2_sym(nind, :))./sqrt(1 + cos(phi2_sym(nind, :)).^2))*mult & Phi_sym(nind, :) <= pi/2 & phi2_sym(nind, :) <= pi/2));
else
    
    % check which are inside the fundamental region
    isInside = phi1_sym <= maxphi1 & phi2_sym <= maxphi2 & Phi_sym <= maxPhi;
end

% take the first one
[~,rowPos] = max(isInside, [], 2);

ind = sub2ind(size(phi1_sym), (1:length(ori))', rowPos);

ori_out = orientation('Euler', [phi1_sym(ind), Phi_sym(ind), phi2_sym(ind)], CS, SS);
end