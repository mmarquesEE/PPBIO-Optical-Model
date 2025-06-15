function dydt = ode_system(~, y, nx, ny, nz, velocity_profile, kon_grid, koff_grid, smax_grid, c0, D_coeff, ru_to_m, dx, dz)
    num_cells = nx * ny * nz;
    c_s = reshape(y(1:num_cells), [nx, ny, nz]);
    s = reshape(y(num_cells + 1:2*num_cells), [nx, ny, nz]);
    %Q = reshape(y(2*num_cells + 1:3*num_cells), [nx, ny, nz]);
    %R = reshape(y(3*num_cells + 1:4*num_cells), [nx, ny, nz]);
    %dcsdt = zeros(nx, ny, nz);
    %dsdt = zeros(nx, ny, nz);

    % Inlet boundary condition (x=1)
    c_s(1,:,:) = c0;
    dcsdt(1,:,:) = 0;

    % Diffusion terms
    d2c_dx2 = zeros(nx, ny, nz);
    d2c_dx2(2:end-1,:,:) = (c_s(3:end,:,:) - 2*c_s(2:end-1,:,:) + c_s(1:end-2,:,:))/ (dx^2);
    
    d2c_dz2 = zeros(nx, ny, nz);
    d2c_dz2(:,:,2:end-1) = (c_s(:,:,3:end) - 2*c_s(:,:,2:end-1) + c_s(:,:,1:end-2))/ (dz^2);
    d2c_dz2(:,:,1) = (c_s(:,:,2) - 2*c_s(:,:,1) + c_s(:,:,1))/ (dz^2);
    d2c_dz2(:,:,end) = (c_s(:,:,end-1) - 2*c_s(:,:,end) + c_s(:,:,end-1))/ (dz^2);
    
    dcsdt = D_coeff * (d2c_dx2 + d2c_dz2);
    
    % Advection
    dcsdt(2:end,:,:) = dcsdt(2:end,:,:) + ...
        bsxfun(@times, velocity_profile, (c_s(1:end-1,:,:) - c_s(2:end,:,:)/ (dx)));
    
    % Adsorption kinetics
    available_sites = max(smax_grid - s, 0);
    dsdt = kon_grid .* c_s .* available_sites - koff_grid .* s;
    dcsdt = dcsdt - (dsdt * ru_to_m);
    
    % Compute dQ/dt and dR/dt
    %dQdt = kon_grid .* c_s + koff_grid;
    %dRdt = c_s .* exp(Q);

    % Combine all derivatives
    %dydt = [dcsdt(:); dsdt(:); dQdt(:); dRdt(:)];
    dydt = [dcsdt(:); dsdt(:)];
end