function [ fx ] = colville(x)
% 
% Colville function 
% Matlab Code by A. Hedar (Sep. 29, 2005).
% The number of variables n = 4.
% 
n = numel(x);
if(~isvector(x))
    error('Input is not a vector');
end

if( n~=4 )
    error('The number of elements in the input is not 4');
end

x1 = x(1);
x2 = x(2);
x3 = x(3);
x4 = x(4);

term1 = 100 * (x1^2-x2)^2;
term2 = (x1-1)^2;
term3 = (x3-1)^2;
term4 = 90 * (x3^2-x4)^2;
term5 = 10.1 * ((x2-1)^2 + (x4-1)^2);
term6 = 19.8*(x2-1)*(x4-1);

fx = term1 + term2 + term3 + term4 + term5 + term6;

% fx  = 100*(x(1)^2-x(2))^2+(x(1)-1)^2+(x(3)-1)^2+90*(x(3)^2-x(4))^2+...
% 10.1*((x(2)-1)^2+(x(4)-1)^2)+19.8*(x(2)^-1)*(x(4)-1);
end
