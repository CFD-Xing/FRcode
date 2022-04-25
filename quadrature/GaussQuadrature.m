function [x, w] = GaussQuadrature(K)
% Define symbolic variable
syms z
% Compute Legendre polynomials
P = LegendrePolynomial(z,K+1);
% Solve Gauss quadrature points [-1, 1]
eqn = P(K) == 0;
x = sort(double(solve(eqn)));
for k=1:K
    z = x(k);
    w(k) = double(subs(2*(1-z^2)/(K+1)^2/P(K+1)^2));
end
end
