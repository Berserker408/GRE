function betahat = lare(x,y,lam,rr,weight)
% Function prla estimates the coefficients of the high dimensional AFT model 
% with the least absolute relative errors (LARE) method. 
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
    rr = 6;
end
it = 0;
delta = ones(p,1);
lambdai = lam;
ynew = log(y);
xnew = x;
betaols = lasso(xnew,ynew,'Lambda',lambdai);
betahat = betaols;
while sum(delta.^2) > toler && it < maxit
    it = it+1;
    beta0 = betahat;
    for j = 1:p
        delta1 = 1;
        it1 = 0;
        while sum(delta1.^2)>toler && it1<maxit1
            betaj = betahat(j);
            derive = sum(weight.*(1-y.^-1.*exp(x*betahat)).*(-y.^-1).*x(:,j).*exp(x*betahat)./(abs(1-y.^-1.*exp(x*betahat))+eps) + ...
                weight.*(1-y.*exp(-x*betahat)).*y.*x(:,j).*exp(-x*betahat)./(abs(1-y.*exp(-x*betahat))+eps)) + ...
                (lambdai*sign(betahat(j))-betahat(j)/rr).*(abs(betahat(j)) <= rr*lambdai);
            hessian = sum(weight.*(1-y.^-1.*exp(x*betahat)).*(-y.^-1).*x(:,j).^2.*exp(x*betahat)./(abs(1-y.^-1.*exp(x*betahat))+eps) + ...
                weight.*(y.^-2).*x(:,j).^2.*exp(2*x*betahat)./(abs(1-y.^-1.*exp(x*betahat))+eps)) + ...
                sum(weight.*(1-y.*exp(-x*betahat)).*(-y).*x(:,j).^2.*exp(-x*betahat)./(abs(1-y.*exp(-x*betahat))+eps) + ...
                weight.*(y.^2).*x(:,j).^2.*exp(-2*x*betahat)./(abs(1-y.*exp(-x*betahat))+eps)) + ...
                (lambdai./abs(betahat(j))-1/rr).*(abs(betahat(j)) <= rr*lambdai);        
            betahat(j) = betahat(j)-derive/hessian;
            it1 = it1+1;
            delta1 = betaj-betahat(j);
        end 
    end
delta = betahat-beta0;
end

