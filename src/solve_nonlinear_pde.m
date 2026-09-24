function [U, info] = solve_nonlinear_pde(epsilon, x, y, bc, options)
%SOLVE_NONLINEAR_PDE Solve u*u_x - u_yy - epsilon*u_xx = f.
% U is Ny-by-Nx; U(j,i) approximates u(x(i),y(j)). x and y include
% endpoints and must be uniformly spaced, increasing vectors.
% bc.left/right: Ny values of u. bc.bottom/top: Nx values of du/dy,
% measured in the POSITIVE y direction on BOTH horizontal edges.
% Dirichlet values take precedence at the four physical corners.
%
% options (all optional):
%   Source            Ny-by-Nx forcing, default zero
%   InitialGuess      Ny-by-Nx warm start, default linear elliptic solve
%   MaxIterations     positive integer, default 200
%   UpdateTolerance   positive relative infinity-norm tolerance, 1e-8
%   ResidualTolerance positive absolute equation tolerance, 1e-8
%   Relaxation        largest trial step in (0,1], default 1
%   MinRelaxation     smallest trial step in (0,Relaxation], default 1/1024
%   Verbose           logical scalar, default false
%
% info.converged must be checked: Picard iteration can stall or reach its
% limit. The central stencil is not guaranteed monotone at high Peclet
% numbers. See docs/method.md for the discretization and its limitations.

    if nargin < 5, options = struct(); end
    validateattributes(epsilon, {'numeric'}, ...
        {'real','finite','scalar','positive'}, mfilename, 'epsilon');
    epsilon = double(epsilon);
    [x, h] = validate_grid(x, 'x');
    [y, p] = validate_grid(y, 'y');
    nx = numel(x); ny = numel(y); n = nx * ny;
    if ~isstruct(bc) || ~isscalar(bc)
        error('pde:Boundary', 'bc must be a scalar structure.');
    end
    names = {'left','right','bottom','top'};
    sizes = [ny, ny, nx, nx];
    for k = 1:numel(names)
        name = names{k};
        if ~isfield(bc, name)
            error('pde:Boundary', 'Missing boundary field: %s.', name);
        end
        validateattributes(bc.(name), {'numeric'}, ...
            {'real','finite','vector','numel',sizes(k)}, mfilename, name);
        value = bc.(name);
        bc.(name) = double(value(:));
    end
    opts = parse_options(options, ny, nx);
    [L, rhs, ids] = diffusion_system(epsilon, h, p, nx, ny, bc, opts.Source);
    if isempty(opts.InitialGuess)
        current = L \ rhs;
    else
        current = reshape(opts.InitialGuess.', [], 1);
    end
    % Exact side values are maintained under relaxation.
    current(1:nx:n) = bc.left;
    current(nx:nx:n) = bc.right;
    if any(~isfinite(current))
        error('pde:InitialSolve', 'The initial state contains nonfinite values.');
    end
    residual = norm((L + convection(current, ids, h, n)) * current - rhs, inf);
    info.converged = residual <= opts.ResidualTolerance;
    info.status = 'maximum_iterations';
    info.history = zeros(0, 4); % iteration, relative update, residual, step
    info.initial_residual = residual;
    info.iterations = 0;
    if info.converged, info.status = 'converged_initial_state'; end

    for iteration = 1:opts.MaxIterations
        if info.converged, break; end
        A = L + convection(current, ids, h, n);
        candidate = A \ rhs;
        if any(~isfinite(candidate))
            info.status = 'nonfinite_linear_solve';
            break;
        end
        linearResidual = norm(A * candidate - rhs, inf);
        linearScale = max(1, norm(A, inf) * norm(candidate, inf) + norm(rhs, inf));
        if linearResidual > 1e-10 * linearScale
            info.status = 'inaccurate_linear_solve';
            break;
        end
        step = opts.Relaxation;
        accepted = false;
        while true
            trial = current + step * (candidate - current);
            trialResidual = norm((L + convection(trial, ids, h, n)) * trial - rhs, inf);
            % Accept a residual decrease, or a value already within tolerance.
            if isfinite(trialResidual) && ...
                    (trialResidual <= residual * (1 - 1e-4 * step) || ...
                     trialResidual <= opts.ResidualTolerance)
                accepted = true;
                break;
            end
            if step <= opts.MinRelaxation, break; end
            step = max(step / 2, opts.MinRelaxation);
        end
        if ~accepted
            info.status = 'line_search_stalled';
            break;
        end
        update = norm(trial - current, inf) / max(1, norm(current, inf));
        current = trial;
        residual = trialResidual;
        info.iterations = iteration;
        info.history(end + 1, :) = [iteration, update, residual, step]; %#ok<AGROW>
        if opts.Verbose
            fprintf('%4d  update %.3e  residual %.3e  step %.4f\n', ...
                iteration, update, residual, step);
        end
        info.converged = update <= opts.UpdateTolerance && ...
                         residual <= opts.ResidualTolerance;
        if info.converged, info.status = 'converged'; end
    end
    U = reshape(current, nx, ny).';
    info.residual_inf = residual;
    info.max_cell_peclet = max(abs(current)) * h / (2 * epsilon);
    info.history_columns = {'iteration','relative_update','residual_inf','relaxation'};
    if info.max_cell_peclet > 1
        warning('pde:UnderresolvedConvection', ...
            'max(|u|)*h/(2*epsilon) = %.3g > 1; refine x or increase epsilon.', ...
            info.max_cell_peclet);
    end
    if ~info.converged
        warning('pde:NotConverged', 'Solver stopped: %s (residual %.3e).', ...
            info.status, info.residual_inf);
    end
end

function [grid, spacing] = validate_grid(grid, name)
    validateattributes(grid, {'numeric'}, ...
        {'real','finite','vector','nonempty'}, mfilename, name);
    grid = double(grid(:));
    if numel(grid) < 3 || any(diff(grid) <= 0)
        error('pde:Grid', '%s must have at least three increasing points.', name);
    end
    spacing = grid(2) - grid(1);
    if max(abs(diff(grid) - spacing)) > 1e-10 * abs(spacing)
        error('pde:Grid', '%s must be uniformly spaced.', name);
    end
end

function opts = parse_options(options, ny, nx)
    if ~isstruct(options) || ~isscalar(options)
        error('pde:Options', 'options must be a scalar structure.');
    end
    opts = struct('Source', zeros(ny,nx), 'InitialGuess', [], ...
        'MaxIterations', 200, 'UpdateTolerance', 1e-8, ...
        'ResidualTolerance', 1e-8, 'Relaxation', 1, ...
        'MinRelaxation', 1/1024, 'Verbose', false);
    fields = fieldnames(options);
    for k = 1:numel(fields)
        name = fields{k};
        if ~isfield(opts, name), error('pde:Options', 'Unknown option: %s.', name); end
        opts.(name) = options.(name);
    end
    validateattributes(opts.Source, {'numeric'}, ...
        {'real','finite','size',[ny,nx]}, mfilename, 'Source');
    if ~isempty(opts.InitialGuess)
        validateattributes(opts.InitialGuess, {'numeric'}, ...
            {'real','finite','size',[ny,nx]}, mfilename, 'InitialGuess');
    end
    validateattributes(opts.MaxIterations, {'numeric'}, ...
        {'scalar','real','finite','integer','positive'}, mfilename, 'MaxIterations');
    positive = {'UpdateTolerance','ResidualTolerance','Relaxation','MinRelaxation'};
    for k = 1:numel(positive)
        validateattributes(opts.(positive{k}), {'numeric'}, ...
            {'scalar','real','finite','positive'}, mfilename, positive{k});
    end
    if opts.Relaxation > 1 || opts.MinRelaxation > opts.Relaxation
        error('pde:Options', 'Require 0 < MinRelaxation <= Relaxation <= 1.');
    end
    validateattributes(opts.Verbose, {'logical'}, {'scalar'}, mfilename, 'Verbose');
    opts.Source = double(opts.Source);
    opts.InitialGuess = double(opts.InitialGuess);
end

function [L, rhs, ids] = diffusion_system(epsilon, h, p, nx, ny, bc, source)
    % Row-major numbering: i + (j-1)*nx. No ghost unknowns are allocated.
    n = nx * ny; rx = epsilon / h^2; ry = 1 / p^2;
    rows = zeros(5*n,1); cols = rows; values = rows; used = 0;
    rhs = reshape(source.', [], 1);
    ids = zeros((nx-2)*ny,1); active = 0;
    for j = 1:ny
        for i = 1:nx
            id = i + (j-1)*nx;
            if i == 1 || i == nx
                used = used + 1;
                rows(used) = id; cols(used) = id; values(used) = 1;
                if i == 1, rhs(id) = bc.left(j); else, rhs(id) = bc.right(j); end
            else
                active = active + 1; ids(active) = id;
                neighbor = [id-1,id,id+1];
                coefficient = [-rx,2*rx+2*ry,-rx];
                if j == 1
                    neighbor(end+1) = id+nx;
                    coefficient(end+1) = -2*ry;
                    rhs(id) = rhs(id) - 2*bc.bottom(i)/p;
                elseif j == ny
                    neighbor(end+1) = id-nx;
                    coefficient(end+1) = -2*ry;
                    rhs(id) = rhs(id) + 2*bc.top(i)/p;
                else
                    neighbor = [neighbor,id-nx,id+nx]; %#ok<AGROW>
                    coefficient = [coefficient,-ry,-ry]; %#ok<AGROW>
                end
                count = numel(neighbor);
                slots = used + (1:count);
                rows(slots) = id; cols(slots) = neighbor; values(slots) = coefficient;
                used = used + count;
            end
        end
    end
    L = sparse(rows(1:used),cols(1:used),values(1:used),n,n);
end

function C = convection(u, ids, h, n)
    coefficient = u(ids) / (2*h);
    C = sparse([ids;ids],[ids-1;ids+1],[-coefficient;coefficient],n,n);
end
