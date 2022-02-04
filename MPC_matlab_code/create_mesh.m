%% Create the static Cartesian grid:

if (strcmp(wave.Fexc_type,'NON_LINEAR_FK'))

    [var] = generate_mesh(X_0, dX, r_cyl, L_cyl);

end

function [var] = generate_mesh(X_0, dX, r_cyl, L_cyl)

    Wx = 1.2 * r_cyl;                         % width of the grid
    Wy = 1.2 * r_cyl;
    Wz = 2.5 * L_cyl;

    Lx = (X_0(1) + Wx) - (X_0(1) - Wx);       % Length of the grid
    Ly = (X_0(2) + Wy) - (X_0(2) - Wy);
    Lz = (X_0(3) + Wz) - (X_0(3) - Wz);

    Nx = floor(Lx/dX(1)); 
    Ny = floor(Ly/dX(2)); 
    Nz = floor(Lz/dX(3)); 

    x_cc = (X_0(1) - Wx) : dX(1) : (X_0(1) + Wx); % cell center coordinates
    y_cc = (X_0(2) - Wy) : dX(2) : (X_0(2) + Wy);
    z_cc = (X_0(3) - Wz) : dX(3) : (X_0(3) + Wz);

    [x, y, z] = ndgrid(x_cc, y_cc, z_cc);         % creating grid points

    X_coord = [ x(:)  y(:)  z(:) ];               % all cell center coordinates

    var = X_coord;

end
