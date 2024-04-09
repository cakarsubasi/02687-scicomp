function result = expspace(base, start, last)
%EXPSPACE Summary of this function goes here
%   Detailed explanation goes here
left = (base*ones(1, last - start));
right = (linspace(start, last, last - start));
result = left .^ right;
end

