function [x, w] = GaussRadauQuadrature(K,s)
% Define symbolic variable
syms z
% Compute Legendre polynomials
P = LegendrePolynomial(z,K);
% Solve Gauss-Radau quadrature points [-1, 1]
% s =  1 -> Left Radau quadrature
% s = -1 -> Right Radau quadrature
eqn = s^K*(P(K) + s*P(K-1))/2 == 0;
x = sort(double(solve(eqn)));
for k=1:K
    z = x(k);
    w(k) = double(subs((1-s*z)/P(K-1)^2/K^2));
end
end
