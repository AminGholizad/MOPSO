function [y] = fitness(x)
y(1) = (sum(cos(3*x),2)).^2;
y(2) = (sum(cos(5*x),2)).^2;
end