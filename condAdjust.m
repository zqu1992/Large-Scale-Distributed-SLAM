function [ y ] = condAdjust( x )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    if (cond(x)>1e4 || sum(sum(x))<1e-10)
        y = eye(length(x));
    else
        y = x;
    end

end

