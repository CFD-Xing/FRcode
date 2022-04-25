function [x, w] = LaguerreGaussQuadrature(K)
% Define symbolic variable
syms z
% Compute Laguerre polynomials
L = LaguerrePolynomial(z,K+1);
% Solve Gauss quadrature points [0, inf[
eqn = L(K) == 0;
x = sort(double(solve(eqn)));
for k=1:K
    z = x(k);
    w(k) = double(subs(z/(K+1)^2/L(K+1)^2));
end
end
