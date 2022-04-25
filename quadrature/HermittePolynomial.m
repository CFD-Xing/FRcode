function [H] = HermittePolynomial(x,K)
% Initialize variable
H = sym(zeros(K+2,1));
% Compute Hermitte polynomials
H(1) = 0; H(2) = 1;
for k=1:K
    H(k+2) = 2*x*H(k+1) - 2*(k-1)*H(k);
end
% Return polynomials
H=H(3:end);
