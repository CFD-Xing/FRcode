function [x, w] = HermitteGaussQuadrature(K)
% Define symbolic variable
syms z
% Compute Hermitte polynomials
H = HermittePolynomial(z,K);
% Solve Gauss quadrature points ]-inf, inf[
eqn = H(K) == 0;
x = sort(double(solve(eqn)));
for k=1:K
    z = x(k);
    w(k) = double(subs(2^(K-1)*factorial(K)*sqrt(pi)/K^2/H(K-1)^2));
end
end
