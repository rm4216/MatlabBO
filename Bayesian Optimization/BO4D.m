clear all
close all

%% Define the grid (sample space)
dimXSpace = 51;
bounds = [-5, 5];
% [X1] = ndgrid(linspace(bounds(1), bounds(2), dimXSpace));
% X = X1(:);
[X1, X2, X3, X4] = ndgrid(linspace(bounds(1), bounds(2), dimXSpace), linspace(bounds(1), bounds(2), dimXSpace),...
    linspace(bounds(1), bounds(2), dimXSpace), linspace(bounds(1), bounds(2), dimXSpace));
X = [X1(:), X2(:), X3(:), X4(:)]; % grid points
[N, d] = size(X);     % # of grid points, dimension
min_x_gap = (bounds(2)-bounds(1))/(dimXSpace-1);

%% Objective function Given
[f, fX, f_min, xs_argmin, f_max ] = ObjectiveFunction( X );

%% Kernel & Hyperparameters for GPs
mu = 0;
meanFunc = {@meanConst};    hyp_f.mean = mu ;
covFunc = {@covSEard};      ell_1 = 1;
                            sf = ceil(f_max/2);   % Prior variance should contain all function values at least within 2*stdDevs when running GP-UCB
                            hyp_f.cov = log([ repmat(ell_1, d, 1); sf ]);
sn = 0.1;
likFunc = @likGauss;        hyp_f.lik = log(sn);

%% Create objective function sampled from GP

%% Compute Hyperparameters from simulation set
noise = normrnd(0, sn, size(fX));
hyp_f = minimize( hyp_f, @gp, -100, @infExact, meanFunc, covFunc, likFunc, X(1:10000:end, :), fX(1:10000:end)+noise(1:10000:end) );
hyp_f.cov(end) = log(sf);       % Keep same prior variance
hyp_f.lik = log(sn);            % Keep same measurement noise
hyp_f.mean = mu;                % Keep same mean
hyp_f.cov(1:end-1) = hyp_f.cov(1:end-1)./5; % reduce regularity assumption, maintain relative relevance

%% Preallocate space for Samples and Observations
maxSamples = dimXSpace^d;
x_sampled = zeros(maxSamples, d);   % Sample all points in the grid
y_observed = zeros(maxSamples, 1);  % Observe all function values in the grid

%% GP-UCB Algorithm
GP_UCB_delta = 0.05;    % Parameter for computing beta_t
delta_f = 0.1;          % Stopping criterion
GP_UCB = 0;             % On/off

%% First Sample
% a=1;
% x_index = floor(a + (N-a).*rand()); % Select initial state uniformly at random from the full grid
x_index = 48;
x_sampled(1, :) = X(x_index, :);
y_observed(1, :) = fX(x_index) + normrnd(0, sn);

%% Make movie in 1 dimension
if(d==1)
    F(maxSamples) = struct('cdata',[],'colormap',[]);
end


    
for k=1:maxSamples
    if(GP_UCB==1)
        
        % Retrain full GP (only lenghtscales) after each iteration
        if(k==1)
            hyp_f = minimize( hyp_f, @gp, -100, @infExact, meanFunc, covFunc, likFunc,...
                x_sampled(1:k, :), y_observed(1:k) );
            hyp_f.cov(end) = log(sf);       % Keep same prior variance
            hyp_f.lik = log(sn);            % Keep same measurement noise
            hyp_f.mean = mu;                % Keep same mean
        end
        
        % Full GP
        [~,~,mean_full,vars_full] = gp(hyp_f, @infExact, meanFunc, covFunc, likFunc,...
            x_sampled(1:k, :), y_observed(1:k), X);
        stdDevs_full = vars_full.^(0.5);
        % beta_t
        beta_t = GP_UCB_beta_t(N, k, GP_UCB_delta);
        % Select new sample
        [~, new_indx] = min(mean_full - beta_t^(0.5).* stdDevs_full); % Try with beta_t reduced by factor 5
        x_star = X(new_indx, :)
        y_x_star = fX(new_indx) + normrnd(0, sn)
        if(d==1)
            frame = Plot1dim(X, fX, mean_full, stdDevs_full, [x_sampled(1:k, :); x_star], [y_observed(1:k); y_x_star]);
            F(k) = getframe(gcf);
        end
        
    else
        if(k==1)
            % Define 1-dim GP
            hyp_1dim.mean = mu;
            hyp_1dim.cov = log([ ell_1; sf ]);
            hyp_1dim.lik = log(sn);
            
            % Define probabilities with softmax of inverse of lengthscales
            pVector = BoltzmanSoftmax(-hyp_f.cov(1:end-1)); % Minus -> inverse of lengthscales used for deriving probability vector over dimensions
        end
        
        explore_dim = DiscreteDistribution(pVector);
        
        % Select lengthscale for 1dim GP
        hyp_1dim.cov(1) = hyp_f.cov(explore_dim);
        
        % Select all samples for 1dim GP
        fixed_components = 1:d ~= repmat(explore_dim, 1, d);
        fixed_values = x_sampled(k, fixed_components);
        
        selection_All = all(X(:, fixed_components) == repmat(fixed_values, N, 1), 2);
        X_1dimAll = X(selection_All, :);
        fX_1dimAll = fX(selection_All);
        
        % Select already sampled x's in that projection
        already_sampled = all(x_sampled(1:k, fixed_components) == repmat(fixed_values, k, 1), 2);
        samples = x_sampled(1:k, :);
        X_1dimSampled = samples(already_sampled, :);
        observations = y_observed(1:k, :);
        Y_observed = observations(already_sampled, :);
        
        % 1dim GP Inference
        [~,~,m_1dim,vars_1dim] = gp(hyp_1dim, @infExact, meanFunc, covFunc, likFunc,...
            X_1dimSampled(:, explore_dim), Y_observed, X_1dimAll(:, explore_dim));
        stdDevs_1dim = vars_1dim.^(0.5);
        
        % Select new sample 1dim GP-UCB
        % beta_t
        beta_t = GP_UCB_beta_t(size(X_1dimAll, 1), length(Y_observed), GP_UCB_delta);
        [~, new_indx] = min(m_1dim - beta_t^(0.5).* stdDevs_1dim);
        x_star = X_1dimAll(new_indx, :)
        y_x_star = fX_1dimAll(new_indx) + normrnd(0, sn)
        
        % Plot 1dim GP and movie frame
        frame = Plot1dim(X_1dimAll(:, explore_dim), fX_1dimAll, m_1dim, stdDevs_1dim, [X_1dimSampled(:, explore_dim); x_star(:, explore_dim)],...
            [Y_observed; y_x_star]);
        %         F(k) = getframe(gcf);
    end
    
    x_sampled(k+1, :) = x_star;
    y_observed(k+1) = y_x_star;
    if(norm(x_star-x_sampled(k, :))<=min_x_gap && abs(y_x_star-y_observed(k))<= 2*sn + delta_f ) % Stopping criterion
        % include && abs(y_x_star-y_best)<= 2*sn
        break;
    end
end

if(d==1)
    frameRate = 2;
    MakeVideoAVI(F(1:k), frameRate);
end

[y_min, y_index] = min(y_observed(1:k));
x_min = x_sampled(y_index, :)
y_min

