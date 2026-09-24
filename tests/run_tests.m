function run_tests()
%RUN_TESTS Analytic invariants, nonzero boundary data, and grid refinement.
% No toolboxes or unit-testing framework required.
    root = fileparts(fileparts(mfilename('fullpath')));
    addpath(fullfile(root,'src'));
    x = linspace(0,1,13); y = linspace(0,1,11);
    [X,Y] = meshgrid(x,y);
    opts = struct('UpdateTolerance',1e-11,'ResidualTolerance',1e-10, ...
        'MaxIterations',300);

    % Zero data exposed the archived ghost-corner division-by-zero defect.
    bc = boundaries(zeros(size(X)),zeros(size(x)),zeros(size(x)));
    [U,info] = solve_nonlinear_pde(0.4,x,y,bc,opts);
    assert(info.converged && max(abs(U(:))) < 1e-12);

    % Nonzero bottom AND top derivative, with non-unit spacing and exact
    % quadratic second differences. This catches Neumann sign/scale errors.
    exact = Y.^2 + 0.3*Y + 0.2*X + 0.1;
    bc = boundaries(exact,0.3*ones(size(x)),2.3*ones(size(x)));
    opts.Source = 0.2*exact - 2;
    [U,info] = solve_nonlinear_pde(0.4,x,y,bc,opts);
    assert(info.converged && max(abs(U(:)-exact(:))) < 1e-8);
    assert(info.residual_inf <= opts.ResidualTolerance);
    limited = opts; limited.MaxIterations = 1;
    [~,stopped] = solve_nonlinear_pde(0.4,x,y,bc,limited);
    assert(~stopped.converged && strcmp(stopped.status,'maximum_iterations'));

    % Smooth manufactured solution: expected second-order global accuracy.
    errors = zeros(1,3);
    for k = 1:3
        n = 2^(k+2)+1;
        x = linspace(0,1,n); y = x; [X,Y] = meshgrid(x,y);
        exact = 0.25*sin(pi*X).*cos(pi*Y) + 0.2*X + 0.1*Y;
        ux = 0.25*pi*cos(pi*X).*cos(pi*Y) + 0.2;
        epsilon = 0.4;
        opts.Source = exact.*ux + (1+epsilon)*0.25*pi^2*sin(pi*X).*cos(pi*Y);
        bc = boundaries(exact,0.1*ones(size(x)),0.1*ones(size(x)));
        [U,info] = solve_nonlinear_pde(epsilon,x,y,bc,opts);
        assert(info.converged);
        errors(k) = max(abs(U(:)-exact(:)));
    end
    orders = log2(errors(1:end-1)./errors(2:end));
    assert(all(orders > 1.7) && all(orders < 2.3));

    % Input failures must be explicit instead of producing a singular solve.
    assert_throws(@() solve_nonlinear_pde(0,x,y,bc,opts));
    bad = bc; bad.left = 0;
    assert_throws(@() solve_nonlinear_pde(0.4,x,y,bad,opts));
    assert_throws(@() solve_nonlinear_pde(0.4,[0,0.2,1],y,bc));
    badOptions = struct('MinRelaxation',0.8,'Relaxation',0.5);
    assert_throws(@() solve_nonlinear_pde(0.4,x,y,bc,badOptions));
    fprintf('All MATLAB checks passed. Refinement orders: %.3f, %.3f\n',orders);
end

function bc = boundaries(U, bottom, top)
    bc = struct('left',U(:,1),'right',U(:,end),'bottom',bottom,'top',top);
end

function assert_throws(action)
    caught = false;
    try, action(); catch, caught = true; end
    assert(caught,'Expected invalid input to raise an error.');
end
