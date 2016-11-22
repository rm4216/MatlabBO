function [ fx ] = camel6( x )
% SIX-HUMP CAMEL FUNCTION
%   Global Minimum
%   fx* = -1.0316
%   x*(1) = (0.0898, -0.7126)
%   x*(2) = (-0.0898, 0.7126)

n = numel(x);
if( n>2 )
    error('The number of elements in the input exceeds 2');
end

x1 = x(1);
x2 = x(2);

term1 = (4-2.1*x1^2+(x1^4)/3) * x1^2;
term2 = x1*x2;
term3 = (-4+4*x2^2) * x2^2;

fx = term1 + term2 + term3;

end

