%% Level-set for tracking a cylinder:

function [psi] = get_level_set_for_cylinder(X1, X_CG, r_cyl, L_cyl)

distance_1 = sqrt( (X1(:,1) - X_CG(1)).^2 + (X1(:,2) - X_CG(2)).^2 ) - r_cyl;  % Level-set calculation for curved section.
distance_2 = X1(:,3) - (X_CG(3) + L_cyl/2);    % Level-set calculation for top surface.
distance_3 = (X_CG(3) - L_cyl/2) - X1(:,3);    % Level-set calculation for bottom surface.

distances = [ distance_1' ; distance_2' ; distance_3' ];

psi = max(distances);

end
