function Z = MCWNNM_ADMM(Y, W, C, rho, mu, K1)
% This routine solves the following weighted nuclear norm optimization problem with column weights,
%
% min_{X, Z} ||W(Y-X)||_F^2 + ||Z||_w,*  s.t.  X = Z
%
% Inputs:
%        Y      -- 3p^2 x M dimensional noisy matrix, D is the data dimension, and N is the number of image patches.
%        NSig -- 3p^2 x 1 dimensional vector of weights
%        Parameters   -- structure of parameters
% Output:
%        Z      -- 3p^2 x M dimensional denoised matrix
% tol = 1e-8;

% Initializing optimization variables
X = zeros(size(Y));
Z = zeros(size(Y));
A = zeros(size(Y));

k = 0;
while k < K1
    k = k + 1;
    
    % update X, fix Z and A
    % min_{X} ||W * Y - W * X||_F^2 + 0.5 * rho * ||X - Z + 1/rho * A||_F^2
    X = diag(1 ./ (W.^2 + 0.5 * rho)) * (diag(W.^2) * Y + 0.5 * rho * Z - 0.5 * A);
    
    % update Z, fix X and A
    % min_{Z} ||Z||_*,w + 0.5 * rho * ||Z - (X + 1/rho * A)||_F^2
    [U, sigma, V] =   svd(full(X + A/rho), 'econ');
    diagSigma = diag(sigma); 
    c1 = diagSigma - eps; 
    c2 = (diagSigma - eps).^2 - 4 * (2/rho*C - eps * diagSigma); 
    c2_positive_index = find(c2 > 0); 
    Z_sigma = max(c1(c2_positive_index) + sqrt(c2(c2_positive_index)), 0)/2; 
    Z = U(:, c2_positive_index) * diag(Z_sigma) * V(:, c2_positive_index)'; 
    
    % update the multiplier A, fix Z and X
    A = A + rho * (X - Z);
    
    rho = min(1e4, mu * rho);
end
