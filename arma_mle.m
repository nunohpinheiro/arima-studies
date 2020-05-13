function results = arma_mle( y, p, q, info )

% ARMA_MLE.M
% results = armae_mle(y,p,q,[info])
% This function computes MLE estimates of the (demeaned) ARMA(p,q) model
%
%  y[t] = phi(1)y[t-1] + ... + phi(p)y[t-p] + e[t] + theta(1)e[t-1] + 
%         ... + theta(q)e[t-q] 
%
% The algorithm computes the exact likelihood function using the Kalman
% filter. It initializes the ARMA coefficients using a method of Hannan and
% Kavalieris.
%
% INPUTS:
% y     : time series of size T with the observations
% p     : autoregressive order of the ARMA model
% q     : moving average order of the ARMA model
% info  : show optimization info =0 no, /=0 yes. If this input is omitted,
%         it is set to 0.
%
% OUTPUT:
% result: structure with the ARMA estimates
%
% REFERENCES: 
% 1) E.J. Hannan and L. Kavalieris. "A Method for Autoregressive-Moving
%       Average Estimation" Biometrika, Vol 71, No2, Aug 1984.
% 2) E.J. Hannan and A.J. McDougall. "Regression Procedures for ARMA
%       Estimation" Journal of the American Statistical Association, Vol
%       83, No 409, June 1988.
% 3) R.H. Jones : Maximum Likelihood Fitting of ARMA Models to
%       Time Series With Missing Observations" Technometrics, Vol 22, No3, 
%       August 1980.
% 
% Written by Constantino Hevia. August 2008

if nargin<4; info = 0; end

% Check data:
[m,n] = size(y); if n>1; T = n; y = y'; else T = m; end
r = max(p,q+1);

% Set optimization options:
if info == 0;
    options = optimset('Largescale','off','Display','off');
else
    options = optimset('Largescale','off','Display','iter');
end

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% INITIALIZE COEFFICIENTS USING Hannan and McDougall (1988) -see below:
coefs = initialize_arma;

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% COMPUTE MLE ESTIMATES USING THE INTERNAL FUNCTION loglikelihood AND THE
% OPTIMIZATION ROUTINE fminunc

coefs1 = fminunc( @log_likelihood, coefs, options );

% Estimate standard deviation of ARMA
[ L, sigmahat, at ] = log_likelihood( coefs1 );

results.ar     = coefs1(1:p);
results.ma     = coefs1(p+1:p+q);
results.sigma  = sigmahat;
results.loglik = -L;

results.residuos = at;                                                      % ATENÇÃO!

    % ---------------------------------------------------------------------
    % INTERNAL FUNCTIONS
  
    function coefs = initialize_arma
        % This function uses Hannan and McDougall (1988) to initialize the
        % coefficients of the ARMA process using a regression approach

        % Step 1: Choose order of autoregression.
        best = 1e+10;
        for h = p : ceil(log(T)^1.5)
            Y = y(h+1:T);
            X = [];
            for j=1:h
                X = [X y(h+1-j : T-j)];
            end
            bet = (X'*X)\X'*Y;
            res = Y - X*bet;
            sigma2 = (res'*res)/length(res);
            BIC = log(sigma2) + h*log(T)/T;
            if BIC < best
                best = BIC; 
                horder=h; 
                residuals = res;
            end
        end
        
        % Step2: Run regression on lagged values and residuals
        nlag = max(p,q);
        Y = y(horder+nlag+1:T);
        X = [];
        for i=1:p 
            X = [X y(horder+nlag+1-i:T-i)];
        end
        for i=1:q
            X = [X residuals(nlag+1-i:T-horder-i)];
        end
        bet = (X'*X)\X'*Y;
        coefs(1:p) = bet(1:p);
        coefs(p+1:p+q) = bet(p+1:p+q);
        coefs = coefs';
    end


    function [ L, sigmahat, at ] = log_likelihood( coefs )
        % function [L sigmahat] = log_likelihood(coefs)
        % This function computes the negative of the log-likelihood of an
        % ARMA(p,q) model:
        %
        %  y[t] = phi(1)y[t-1] + ... + phi(p)y[t-p] + ... 
        %         + e[t] + theta(1)e[t-1] + ... + theta(q)e[t-q] 
        %
        % where phi    = coefs(1:p)
        %       theta  = coefs(p+1:p+q)
        
        phi   = coefs(1:p); 
        theta = coefs(p+1:p+q);
        
        % Build matrices of state space representation:
        %
        %      x[t+1] = A x[t] + R eps[t+1]
        %      y[t]   = Z'x[t]
        %
        
        A = zeros(r,r); R = zeros(r,1); Z = zeros(r,1);
        A(1:p,1)=phi; A(1:r-1,2:r)=eye(r-1);
        R(1) = 1; R(2:q+1) = theta;
        Z(1) = 1; Q = real( R*R' );                                         % METI REAL()
        
        % Initialize state and covariance matrix
        xhat_tt  = zeros(r,1);
        Sigma_tt = reshape( (eye(r^2)-kron(A,A) ) \ Q(:), r, r );
        
        logsum = 0;
        sumsq  = 0;
        
        at = 0;                                                             % ATENÇÃO!

        for i=0:T-1     % Start Kalman Filter Recursion
            xhat_t1t  = real( A*xhat_tt );                                  % METI REAL()
            Sigma_t1t = real( A*Sigma_tt*A' + Q );                          % METI REAL() 
            yhat_t1t  = Z'*xhat_t1t;
            omega     = Z'*Sigma_t1t*Z;
            delta_t1  = Sigma_t1t*Z/omega;
            innov     = y(i+1) - yhat_t1t;
            xhat_t1t1 = xhat_t1t + delta_t1*innov;
            Sigma_t1t1= Sigma_t1t - delta_t1*Z'*Sigma_t1t;
            
            % Add likelihood terms:
            logsum = logsum + log(omega);
            sumsq  = sumsq  + innov^2/omega;
            
            at = [at sqrt(innov^2/omega)];                                  % ATENÇÃO!
            
            % Update estimates;
            xhat_tt  = xhat_t1t1; 
            Sigma_tt = Sigma_t1t1;
        end
        
        L = real( logsum + T*log(sumsq) );                                  % METI REAL()
        sigmahat = real( sqrt(sumsq/T) );                                   % METI REAL()
    
    end % function log_likelihood


end %function arma_mle
    

