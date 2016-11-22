function [ beta_t ] = GP_UCB_beta_t( D, t, delta )

% Computation of beta_t, Theorem1 GP-UCB Bounds for GP prior for finite D

beta_t = 2 * log( (D * t^2 * pi^2) / (6*delta) );

end