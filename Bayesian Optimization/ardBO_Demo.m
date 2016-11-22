clear all
close all

%% Define the grid (sample space)
dimXSpace = 101;
n=2;
[X1, X2] = ndgrid(linspace(-2, 2, dimXSpace), linspace(-1, 1, dimXSpace));
X = [X1(:), X2(:)]; % grid points
N = size(X, 1);     % # of grid points

%% Objective function
f = @camel6;

%% Plot objective function
fX = zeros(N, 1);
for i = 1 : N
    fX(i) = f(X(i, :));
end
fX1X2 = reshape(fX, [dimXSpace, dimXSpace]);
figure();
mesh(X1, X2, fX1X2);

%% f* x*
f_star = min(fX);
indices = find( fX == f_star );
xs_star = X(indices, :);

%% Kernel & Hyperparameters for GP-Learning
k = 0;
meanFunc = {@meanConst};    hyp_f.mean = k ;
covFunc = {@covSEard};      ell_1 = 1; ell_2 = 1; sf = 5; hyp_f.cov = log([ ell_1; ell_2; sf ]);
sn = 0.1;
likFunc = @likGauss; hyp_f.lik = log(sn);

hyp_1dim.mean = k;
hyp_1dim.cov = log([ ell_1; sf ]);
hyp_1dim.lik = log(sn);

%% Initialization
maxSamples = dimXSpace^2;
x_sampled = zeros(maxSamples, 2);
y_observed = zeros(maxSamples, 1);

%% ardBO algorithm
a=1;
x_index = floor(a + (N-a).*rand());

x_sampled(1, :) = X(x_index, :);
y_observed(1, :) = f(x_sampled(1, :)) + normrnd(0, sn);

for k = 1:maxSamples
    if(k==1)
        hyp_f = minimize( hyp_f, @gp, -100, @infExact, meanFunc, covFunc, likFunc,...
            x_sampled(1:k, :), y_observed(1:k) );
    end
    [ ~ , explore_indx ] = min(hyp_f.cov(1:end-1));
    
    % Parameters for 1-dim GP
    hyp_1dim.cov(1) = hyp_f.cov(explore_indx);
%     hyp_1dim.cov = hyp_f.cov([explore_indx, end]);
%     hyp_1dim.lik = hyp_f.lik;
    
    % Select samples for 1dimGP
    fixed_components = 1:n ~= [explore_indx, explore_indx];
    fixed_value = x_sampled(k, fixed_components);
    selection_All = X(:, fixed_components)==fixed_value*ones(N, 1); % may need to specify à of fixed components instead of 1
    X_1dimAll = X(selection_All, :);
    
    % Select already sampled x's in that projection
    already_visited = x_sampled(1:k, fixed_components) == fixed_value;
    samples = x_sampled(1:k, :);
    X_1dimVisited = samples(already_visited, :);
    observations = y_observed(1:k, :);
    Y_visited = observations(already_visited, :);
    
    
    [~,~,m_1dim,vars_1dim] = gp(hyp_1dim, @infExact, meanFunc, covFunc, likFunc,...
                                X_1dimVisited(:, explore_indx), Y_visited, X_1dimAll(:, explore_indx));
    
    stdDevs_1dim = sqrt(vars_1dim);
    figure(44); clf
    plot(X_1dimAll(:, explore_indx), m_1dim); hold on
    plot(X_1dimAll(:, explore_indx), m_1dim+2*stdDevs_1dim)
    hold on
    plot(X_1dimAll(:, explore_indx), m_1dim)
    plot(X_1dimAll(:, explore_indx), m_1dim-2*stdDevs_1dim)
    
    % Select new sample
    t = length(Y_visited);
    delta = 0.05;
    pi_t = (pi^2) * (t^2) / 6;
    beta_t = 2 * log(size(X_1dimAll, 1)*pi_t/delta);
    
    [~, new_indx] = min(m_1dim - beta_t.* stdDevs_1dim);
%     x_sampled(k+1, fixed_components) = fixed_value;
%     x_sampled(k+1, ~fixed_components)= new_value;
    x_sampled(k+1, :) = X_1dimAll(new_indx, :);
    y_observed(k+1) = f(x_sampled(k+1, :)) + normrnd(0, sn);
    hyp_current_dim = minimize( hyp_1dim, @gp, -100, @infExact, meanFunc, covFunc, likFunc,...
                        [X_1dimVisited(:, explore_indx);x_sampled(k+1, explore_indx)], [Y_visited; y_observed(k+1)] );
    hyp_1dim.cov(1) = hyp_current_dim.cov(1);
    hyp_f.cov(explore_indx) = hyp_1dim.cov(1);
    
    selection_Alternative = X(:, ~fixed_components)== x_sampled(k+1, ~fixed_components);
    X_1dimAlternative = X(selection_Alternative, :);
    
end

