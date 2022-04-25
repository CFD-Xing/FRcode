%% Clean-up
clc
clear all
close all

%% Gauss Quadrature
fprintf('------------------------------\n');
fprintf('GAUSS QUADRATURE\n')
fprintf('------------------------------\n');
for K=2:5
    fprintf('ORDER %i\n', K)
    [x, w] = GaussQuadrature(K);
    for i=1:K
        fprintf('  x%i = %1.6e ', i, x(i))
        fprintf('  w%i = %1.6e\n', i, w(i))
    end
    fprintf('\n')
end
fprintf('\n')

%% Gauss-Radau Quadrature (Left)
fprintf('------------------------------\n');
fprintf('GAUSS-RADAU QUADRATURE (LEFT)\n');
fprintf('------------------------------\n');
for K=2:5
    fprintf('ORDER %i\n', K)
    [x, w] = GaussRadauQuadrature(K,1);
    for i=1:K
        fprintf('  x%i = %1.6e ', i, x(i))
        fprintf('  w%i = %1.6e\n', i, w(i))
    end
    fprintf('\n')
end
fprintf('\n')

%% Gauss-Radau Quadrature (Right)
fprintf('------------------------------\n');
fprintf('GAUSS-RADAU QUADRATURE (RIGHT)\n');
fprintf('------------------------------\n');
for K=2:5
    fprintf('ORDER %i\n', K)
    [x, w] = GaussRadauQuadrature(K,-1);
    for i=1:K
        fprintf('  x%i = %1.6e ', i, x(i))
        fprintf('  w%i = %1.6e\n', i, w(i))
    end
    fprintf('\n')
end
fprintf('\n')

%% Gauss-Lobatto Quadrature
fprintf('------------------------------\n');
fprintf('GAUSS-LOBATTO\n');
fprintf('------------------------------\n');
for K=3:5
    fprintf('ORDER %i\n', K)
    [x, w] = GaussLobattoQuadrature(K);
    for i=1:K
        fprintf('  x%i = %1.6e ', i, x(i))
        fprintf('  w%i = %1.6e\n', i, w(i))
    end
    fprintf('\n')
end
fprintf('\n')

%% Hermitte-Gauss Quadrature
fprintf('------------------------------\n');
fprintf('HERMITTE-GAUSS QUADRATURE\n')
fprintf('------------------------------\n');
for K=2:5
    fprintf('ORDER %i\n', K)
    [x, w] = HermitteGaussQuadrature(K);
    for i=1:K
        fprintf('  x%i = %1.6e ', i, x(i))
        fprintf('  w%i = %1.6e\n', i, w(i))
    end
    fprintf('\n')
end
fprintf('\n')

%% Laguerre-Gauss Quadrature
fprintf('------------------------------\n');
fprintf('LAGUERRE-GAUSS QUADRATURE\n')
fprintf('------------------------------\n');
for K=2:5
    fprintf('ORDER %i\n', K)
    [x, w] = LaguerreGaussQuadrature(K);
    for i=1:K
        fprintf('  x%i = %1.6e ', i, x(i))
        fprintf('  w%i = %1.6e\n', i, w(i))
    end
    fprintf('\n')
end
fprintf('\n')

%% Chebyshev-Gauss Quadrature
fprintf('------------------------------\n');
fprintf('CHEBYSHEV-GAUSS QUADRATURE\n')
fprintf('------------------------------\n');
for K=2:5
    fprintf('ORDER %i\n', K)
    [x, w] = ChebyshevGaussQuadrature(K);
    for i=1:K
        fprintf('  x%i = %1.6e ', i, x(i))
        fprintf('  w%i = %1.6e\n', i, w(i))
    end
    fprintf('\n')
end
fprintf('\n')
