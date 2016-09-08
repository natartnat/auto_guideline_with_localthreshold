function [ pout ] = sortedPoint( p )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

px = p(:,1);
py = p(:,2);

cx = mean(px);
cy = mean(py);

a = atan2(py-cy,px-cx);

[~,ord] = sort(a);

px = px(ord);
py = py(ord);

pout = [px py];

end