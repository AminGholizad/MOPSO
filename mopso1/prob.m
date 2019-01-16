function [y,c] = prob(x)
y(1) = (sum(cos(3*x),2)).^2;
y(2) = (sum(cos(5*x),2)).^2;
c = max(x(1) - x(2) + eps, 0) + max(x(2) - x(3) + eps,0) + max((sum(cos(x), 2) - 0.8)^4 - 1e-5 + eps, 0);
end