% DEMO_INTERNAL_LAYER A modest-grid, self-contained portfolio demonstration.
% Run from any working directory. MATLAB base functionality only.
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root, 'src'));
x = linspace(0,6,121); y = linspace(0,5,31);
bc.left = ones(numel(y),1); bc.right = -ones(numel(y),1);
bc.bottom = zeros(1,numel(x)); bc.top = zeros(1,numel(x));
epsilon_sequence = [0.2, 0.1, 0.05];
opts = struct('Verbose',true,'MaxIterations',300, ...
    'UpdateTolerance',1e-8,'ResidualTolerance',1e-8);
profiles = zeros(numel(epsilon_sequence),numel(x));
fprintf(' epsilon    iterations    residual_inf    cell_Peclet\n');
for k = 1:numel(epsilon_sequence)
    [U, info] = solve_nonlinear_pde(epsilon_sequence(k),x,y,bc,opts);
    if ~info.converged
        error('demo:NotConverged', 'epsilon=%g: %s',epsilon_sequence(k),info.status);
    end
    profiles(k,:) = U(ceil(numel(y)/2),:);
    fprintf('%8.3g    %5d         %.3e       %.3f\n', ...
        epsilon_sequence(k),info.iterations,info.residual_inf,info.max_cell_peclet);
    opts.InitialGuess = U;
end
figure('Color','w');
subplot(1,2,1);
surf(x,y,U,'EdgeColor','none'); view(40,30); colorbar;
xlabel('x'); ylabel('y'); zlabel('u'); title('Regularized steady state');
subplot(1,2,2);
plot(x,profiles,'LineWidth',1.5); grid on;
xlabel('x'); ylabel('u(x, mid-y)'); title('Viscosity continuation');
legend(arrayfun(@(v) sprintf('epsilon = %.2g',v),epsilon_sequence, ...
    'UniformOutput',false),'Location','southwest');
% For export, save the figure manually or use exportgraphics in recent MATLAB.
