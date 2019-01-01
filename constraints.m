function y = constraints(x)
y(1) = x(1)<x(2);
y(2) = x(2)<x(3);
y(3) = (sum(cos(x),2) - 0.8)^4 < 1e-5;
end

