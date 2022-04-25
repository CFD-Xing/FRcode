function [T] = ChebyshevPolynomial(x,K)
% Initialize variable
T = sym(zeros(K+1,1));
% Compute Chebyshev polynomials
T(1) = 1; T(2) = x;
for k=1:K-1
    T(k+2) = 2*x*T(k+1) - T(k);
end
% Return polynomials
T=T(2:end);
