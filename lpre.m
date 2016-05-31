function betahat = lpre(x,y,lam,rr,weight)
% Function lpre estimates the coefficients of the high dimensional AFT model 
% with the least product relative errors (LPRE) method. 
% We use the MCP penalty in this function, we also can change it with Lasso 
% or SCAD.
% Input:
%	x: n x p covariates.
% 	y: n dimensional response vector.
%	lam, rr: tuning parameters.
%	weight: the weight allocated to the samples.  n dimensional vector. 
%           Default as 1/n. n is the number of samples.         
% Ouput:
%   betahat: estimated coefficient. p dimensional vector.
% Yangguang Zang <yangguang.zang@gmail.com>
% $Revision: 1.0.0 $  $Date: 2016/05/03 $


maxit = 50;
maxit1 = 50;
toler = 1E-6;
[n, p] = size(x);
if nargin < 5
    weight = ones(n,1)/n;
end
if nargin < 4
    rr  =  6;
end
it = 0;
delta = ones(p,1);
ynew = log(y);
xnew = x;
betaols = lasso(xnew,ynew,'Lambda',lam);
betahat = betaols;
while sum(delta.^2) > toler && it < maxit
    it = it+1;
    beta0 = betahat;
    for j = 1:p
        delta1 = 1;
        it1 = 0;
        while sum(delta1.^2) > toler && it1 < maxit1
            betaj = betahat(j);
            derive = sum(weight.*exp(x*betahat).*y.^(-1).*x(:,j) - ...
                weight.*y.*exp(-x*betahat).*x(:,j)) + ...
                (lam*sign(betaj)-betaj/rr).*(abs(betaj) <= rr*lam);
            hessian = sum(weight.*exp(2*x*betahat).*y.^(-1).*x(:,j).^2 + ...
                weight.*y.*exp(-2*x*betahat).*x(:,j).^2) + ...
                ((lam./(abs(betaj))-1/rr).*(abs(betaj) <= rr*lam));        
            betahat(j) = betahat(j) - derive/hessian;
            it1 = it1 + 1;
            delta1 = betaj - betahat(j);
        end 
    end
    delta = betahat-beta0;
end
    
    