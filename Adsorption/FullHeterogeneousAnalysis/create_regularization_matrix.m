function L = create_regularization_matrix(nx, ny)
    % Creates a sparse matrix L that approximates the gradient operator 
    % for a 2D grid of size [nx, ny] that has been vectorized.
    
    % 1D derivative matrix for ny points (Y-direction, along columns)
    Dy = spdiags([-ones(ny,1), ones(ny,1)], [0, 1], ny, ny);
    % Set the last entry to zero to prevent wrapping
    Dy(ny,ny) = 0; 
    
    % 1D derivative matrix for nx points (X-direction, along rows)
    Dx = spdiags([-ones(nx,1), ones(nx,1)], [0, 1], nx, nx);
    Dx(nx,nx) = 0;
    
    % Identity matrices for Kronecker product
    Ix = speye(nx);
    Iy = speye(ny);
    
    % Create 2D operators using the Kronecker product
    % Ly calculates differences between vertical neighbors
    % Lx calculates differences between horizontal neighbors
    Ly = kron(Ix, Dy);
    Lx = kron(Dx, Iy);
    
    % Combine to form the full gradient operator for a single parameter map
    L_single_param = [Lx; Ly];
    
    % Create a block diagonal matrix for all 3 parameter types (kon, koff, smax)
    L = blkdiag(L_single_param, L_single_param, L_single_param);
end