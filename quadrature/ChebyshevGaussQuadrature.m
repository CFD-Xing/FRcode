function [x, w] = ChebyshevGaussQuadrature(K)
% Define symbolic variable
syms z
% Compute Chebyshev polynomials
T = ChebyshevPolynomial(z,K);
% Solve Gauss quadrature points [-1, 1]
eqn = T(K) == 0;
x = sort(double(solve(eqn)));
for k=1:K
    z = x(k);
    w(k) = pi/K;
end
end
