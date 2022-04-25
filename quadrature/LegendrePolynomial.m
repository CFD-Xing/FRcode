function [P] = LegendrePolynomial(x,K)
% Initialize variable
P = sym(zeros(K+2,1));
% Compute Legendre polynomials
P(1) = 0; P(2) = 1;
for k=1:K
    P(k+2) = (2*k-1)/k*x*P(k+1) - (k-1)/k*P(k);
end
% Return polynomials
P=P(3:end);