function [ f, fX, f_min, xs_argmin, f_max ] = ObjectiveFunction( X )

[N, d] = size(X);
f = @levy;

%% Evaluate objective function
fX = EvaluateAll(X, f, N, d);
if(d==1)
    figure();
    plot(X, fX);
end

%% f* x*
f_min = min(fX);
indices =  fX == f_min ;
xs_argmin = X(indices, :);
f_max = max(fX); % Necessary for defining prior variance "sf" when running GP-UCB

end

