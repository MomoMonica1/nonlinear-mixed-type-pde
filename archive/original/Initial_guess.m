function U_initial = Initial_guess(epsilon, h, p, Nx, Ny, bottom, top, left, right)
%INITIAL_GUESS function to compute the solution of Uyy+εUxx=0 as 
% initial guess
%
% Input arguments:
%       epsilon - Parameter of the PDE
%       h - Grid spacing in x direction
%       p - Grid spacing in y direction 
%       Nx - number of gird points in x direction
%       Ny - number of grid points in y direction
%       bottom - bottom boundary condition
%       top - top boundary condition
%       left - left boundary condition
%       right - right boundary condition
%
% Output argument:
%       U_initial - solution vector of the initial guess

%Set up the main matrix
    rx = epsilon/(h^2);
    ry = 1/(p^2);
    e = ones(Nx*Ny,1);
    maindiag = (2*rx+2*ry)*e;
    xdiag = -rx*e;
    ydiag = -ry*e;
    
    A = spdiags([ydiag, xdiag, maindiag, xdiag, ydiag], ...
                [-Nx, -1, 0, 1, Nx], Nx*Ny, Nx*Ny);

    F = zeros(Nx*Ny,1);
    
    %Create the global-index matrix
    [index_x, index_y] = meshgrid(1:Nx, 1:Ny);
    id_mat = (index_y-1)*Nx + index_x;
    
    %collect the global-indices of boundary points 
    leftbd = id_mat(:,1);
    rightbd = id_mat(:,end);
    bottombd = id_mat(1,:);
    topbd = id_mat(end,:);
    
   % Insert boundary conditions into the main matrix
   % Lower boundary
    for k = 2:(length(bottombd)-1)
        idx = bottombd(k);
        A(idx,:) = 0;
        A(idx,idx+Nx) = 1;
        A(idx,idx) = -1;
        F(idx) = bottom(k);
    end

    % Upper boundary
    for k = 2:(length(topbd)-1)
        idx = topbd(k);
        A(idx,:) = 0;
        A(idx,idx) = -1;
        A(idx,idx-Nx) = 1;
        F(idx) = top(k);
    end
    
    % Left boundary
    for k = 1:length(leftbd)
        idx = leftbd(k);
        A(idx,:) = 0;
        A(idx,idx) = 1;
        F(idx) = left(k);
    end

    % Right boundary
    for k = 1:length(rightbd)
        idx = rightbd(k);
        A(idx,:) = 0;
        A(idx,idx) = 1;
        F(idx) = right(k);
    end
    
    % Solve the linear system
    U_initial = A\F;
end