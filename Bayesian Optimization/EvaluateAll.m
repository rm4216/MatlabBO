function [ fX ] = EvaluateAll( X, f, NN, dd)
% Evaluate f at all points in a domain X
% f:R^dim->R
% N data points
% Return a vector of N function values

[N, dim] = size(X);

if(all([NN, dd]==[N, dim]))
    fX = zeros(N, 1);
    for i=1:N
        fX(i) = f(X(i, :));
    end
else
    error('Unexpected X size');
end

end

