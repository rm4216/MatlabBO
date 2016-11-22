function [ index ] = DiscreteDistribution( pVector )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if(~isvector(pVector))
    error('Input probabilities should be a vector');
end

r = rand;

index = sum(~(r<cumsum(pVector))) + 1;

end

