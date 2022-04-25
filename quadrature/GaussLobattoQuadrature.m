function [x, w] = GaussLobattoQuadrature(K)
% Define symbolic variable
syms z
% Compute Legendre polynomials
P = LegendrePolynomial(z,K);
% Solve Gauss-Lobatto quadrature points [-1, 1]
eqn = z*P(K-1) - P(K-2) == 0;
x = sort(double(solve(eqn)));
for k=1:K
    z = x(k);
    w(k) = double(subs(2/P(K-1)^2/K/(K-1)));
end
end
