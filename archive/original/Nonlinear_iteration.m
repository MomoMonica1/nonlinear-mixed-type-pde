function U_final = Nonlinear_iteration(U_initial, epsilon, h, p, Nx, Ny, ...
                                      bottom, top, left, right, max_iterations, error_tolerance)
% NONLINEAR_ITERATION Solve 2D PDE U*Ux - Uyy - eps*Uxx = 0 with FD + iteration
% --- Inputs as in your original description ---
% --- Outputs: U_final (column vect
% or of interior solution, size Nx*Ny) ---

    % coefficients
    rx = epsilon/(h^2);
    ry = 1/(p^2);

    Ny_ghost = Ny + 2;                 % total rows including ghost rows
    total_points = Nx * Ny_ghost;      % total unknowns in the ghosted grid

    % index matrix: linear indices arranged row-by-row (y rows, x columns)
    [x_index, y_index] = meshgrid(1:Nx, 1:Ny_ghost);
    index_matrix = (y_index-1)*Nx + x_index;

    % initialize U_current: ghosted array of size (Ny_ghost x Nx)
    U_current = zeros(Ny_ghost, Nx);
    % U_initial is interior with size Nx*Ny, reshape into (Ny x Nx)
    U_initial_mat = reshape(U_initial, Nx, Ny)';  % Ny x Nx
    % place interior into rows 2..Ny_ghost-1
    U_current(2:Ny_ghost-1, :) = U_initial_mat;

    % precompute boundaries indices (column-major linear indices)
    left_boundary  = index_matrix(:,1);       % all rows, first column
    right_boundary = index_matrix(:,end);     % all rows, last  column
    bottom_boundary = index_matrix(2,:);      % interior bottom row (row 2)
    bottom_ghost    = index_matrix(1,:);      % ghost bottom   (row 1)
    top_boundary    = index_matrix(end-1,:);  % interior top    (row Ny_ghost-1)
    top_ghost       = index_matrix(end,:);    % ghost top       (row Ny_ghost)

    % main iteration loop
    for iteration = 1:max_iterations

        % Build linear operator for -eps*Uxx - Uyy (5-diagonal)
        e = ones(total_points,1);
        main_diagonal = (2*rx + 2*ry) * e;
        x_diagonal = -rx * e;
        y_diagonal = -ry * e;

        A = spdiags([y_diagonal, x_diagonal, main_diagonal, x_diagonal, y_diagonal], ...
                    [-Nx, -1, 0, 1, Nx], total_points, total_points);
        F = zeros(total_points,1);

        % ---- Dirichlet on left and right boundaries (apply to all interior rows) ----
        % left and right are assumed length Ny (interior y points). We map them
        % to rows 2..Ny_ghost-1.
        for row = 2:Ny_ghost-1
            idxL = index_matrix(row, 1);
            idxR = index_matrix(row, Nx);

            % left boundary (Dirichlet)
            A(idxL, :) = 0;
            A(idxL, idxL) = 1;
            F(idxL) = left(row-1);    % left(1..Ny) map to row 2..Ny_ghost-1

            % right boundary (Dirichlet)
            A(idxR, :) = 0;
            A(idxR, idxR) = 1;
            F(idxR) = right(row-1);
        end

        % ---- Bottom ghost boundary (Neumann) ----
        % Using central difference: (U_2 - U_0) / (2p) = g_bottom  => U_0 = U_2 - 2p*g
        for col = 1:Nx
            id_ghost = bottom_ghost(col);     % row 1
            id_row2  = bottom_boundary(col);  % row 2 (interior)

            A(id_ghost, :) = 0;
            A(id_ghost, id_ghost) = 1;       % U_0
            A(id_ghost, id_row2)  = -1;      % -U_2
            F(id_ghost) = -2 * p * bottom(col);  % -2p*g_bottom(col)
        end

        % ---- Top ghost boundary (Neumann) ----
        % Using central difference at top: (U_N+1 - U_N-1)/(2p) = g_top
        % => U_ghost = U_{row-2} + 2p*g_top  (row-2 is two rows below ghost)
        for col = 1:Nx
            id_ghost = top_ghost(col);            % row Ny_ghost
            id_row_minus2 = index_matrix(end-2, col); % row Ny_ghost-2

            A(id_ghost, :) = 0;
            A(id_ghost, id_ghost) = 1;
            A(id_ghost, id_row_minus2) = -1;
            F(id_ghost) = 2 * p * top(col);
        end

        % ---- Corner ghost / boundary assignments ----
        % (Left-top, right-top, left-bottom, right-bottom) -- keep as placeholders
        % Replace these with correct physical BCs if needed.
        % bottom-left ghost (row 1, col 1)
        A(index_matrix(1,1), :) = 0;
        A(index_matrix(1,1), index_matrix(1,1)) = left(1);
        F(index_matrix(1,1)) = 1;  % placeholder

        % bottom-right ghost (row 1, col Nx)
        A(index_matrix(1,Nx), :) = 0;
        A(index_matrix(1,Nx), index_matrix(1,Nx)) = right(1);
        F(index_matrix(1,Nx)) = -1; % placeholder

        % top-left ghost (row Ny_ghost, col 1)
        A(index_matrix(end,1), :) = 0;
        A(index_matrix(end,1), index_matrix(end,1)) = left(end);
        F(index_matrix(end,1)) = 1; % placeholder

        % top-right ghost (row Ny_ghost, col Nx)
        A(index_matrix(end,Nx), :) = 0;
        A(index_matrix(end,Nx), index_matrix(end,Nx)) = right(end);
        F(index_matrix(end,Nx)) = -1; % placeholder

        % ---- Add nonlinear term U * U_x (using central difference in x) ----
        % We only modify interior points where both left and right neighbors exist:
        % columns 2 .. Nx-1, rows 2 .. Ny_ghost-1
        U_current_vec = reshape(U_current.', [], 1); % column vector size total_points
        U_mat = reshape(U_current_vec, Nx, Ny_ghost)'; % Ny_ghost x Nx, consistent layout

        for row = 2:Ny_ghost-1
            for col = 2:Nx-1
                idx = index_matrix(row, col);
                Uval = U_mat(row, col);  % U at the grid point
                % central FD: U_x ~ (U_{i+1} - U_{i-1})/(2*h)
                % contribution to matrix (move U*Ux terms to LHS): modify off-diagonals
                A(idx, idx-1) = A(idx, idx-1) - Uval/(2*h);
                A(idx, idx+1) = A(idx, idx+1) + Uval/(2*h);
            end
        end

        % ---- Solve linear system ----
        U_new = A \ F;  % total_points x 1

        % ---- compute relative error over all interior points (rows 2..Ny_ghost-1, all cols) ----
        interior_indices = index_matrix(2:Ny_ghost-1, 1:Nx); % Ny x Nx
        interior_id_vec = interior_indices';                  % transpose to Nx x Ny
        interior_id = interior_id_vec(:);                     % column vector of length Nx*Ny

        % Ensure U_current_vec is the vector form of current iterate
        % (U_current_vec was created earlier)
        rel_num = norm(U_new(interior_id) - U_current_vec(interior_id));
        rel_den = norm(U_current_vec(interior_id)) + 1e-8;
        relative_error = rel_num / rel_den;

        fprintf('Iteration %d, relative error = %.4e\n', iteration, relative_error);

        if relative_error < error_tolerance
            fprintf('Converged at iteration %d, relative error = %.4e\n', iteration, relative_error);
            U_current = U_new; % store the converged vector
            break;
        end

        % update iterate
        U_current = U_new; % vector form is fine; in next iteration we reshape as needed
    end

    % Final arrangement: reshape U_current (vector) back to ghosted matrix, extract interior
    U_final_mat = reshape(U_current, Nx, Ny_ghost)';   % Ny_ghost x Nx
    interior_U = U_final_mat(2:end-1, :);               % Ny x Nx
    U_final = interior_U';                              % Nx x Ny
    U_final = U_final(:);                               % column vector (Nx*Ny x 1)

    if iteration == max_iterations
        fprintf('Reached maximum iterations, final relative error = %.4e\n\n', relative_error);
    end
end