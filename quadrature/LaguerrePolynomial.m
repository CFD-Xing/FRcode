function [L] = LaguerrePolynomial(x,K)
% Initialize variable
L = sym(zeros(K+2,1));
% Compute Laguerre polynomials
L(1) = 0; L(2) = 1;
for k=1:K
    L(k+2) = (2*k-1-x)/k*L(k+1) - (k-1)/k*L(k);
end
% Return polynomials
L=L(3:end);
