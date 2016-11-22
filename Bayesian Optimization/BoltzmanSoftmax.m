function [ pVector ] = BoltzmanSoftmax( vVector )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
if(~isvector(vVector))
    error('Input values should be a vector');
end

pVector = exp(vVector)./sum(exp(vVector));

end

