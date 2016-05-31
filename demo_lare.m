% This Demo shows how to use functions lare/lare to estimate the coefficients  
% of the high dimensional AFT model. We also use the 5-fold cross validation
% method to select the best lambda.
% Yangguang Zang <yangguang.zang@gmail.com>
% $Revision: 1.0.0 $  $Date: 2016/05/03 $
tic,
clear
n = 100;
p = 200;
maxit = 20;
maxit1 = 20;
toler = 1E-4;
nLambda = 20;
rr = 6;
lambdaRatio = 1E-5;
%%%%%%%%%%%%%%%%%%%%%%%%
w = randn(n,p);
wtest = randn(n,p);
x = w;
%%%%%%%%%%%%%%%%%%%%%%%%%
beta = zeros(p,1);
beta(1:10) = 0.4+0.8*rand(10,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%
y = exp(x*beta).*exp(randn(n,1));
beta0 = zeros(p,1);
for i = 1:100
    lam = i;
    betahat = lare(x,y,lam,6);
    if sum(abs(betahat) <= 0.001) == p
        break;
    end
end
lambdaMax = i;
lambdaMin = lambdaMax * lambdaRatio;
loghi = log(lambdaMax);
loglo = log(lambdaMin);
logrange = loghi - loglo;
interval = -logrange/(nLambda-1);
lambda = exp(loghi:interval:loglo)';
betaset = zeros(p,nLambda);
bicset = zeros(nLambda,1);
for i = 1:nLambda
    i
    for ttt = 1:5
        xtest = x(floor((ttt-1)*n/5)+1:floor(ttt*n/5),:);
        ytest = y(floor((ttt-1)*n/5)+1:floor(ttt*n/5),:);
        xtrain = x(setdiff(1:n,floor((ttt-1)*n/5)+1:floor(ttt*n/5)),:);
        ytrain = y(setdiff(1:n,floor((ttt-1)*n/5)+1:floor(ttt*n/5)),:);        
    	beta0 = lare(xtrain,ytrain,lambda(i),6);
        bicset(i,1) = bicset(i,1)+sum((log(ytest)-xtest*beta0).^2);
    end
end
bicbesti = find(bicset == min(bicset(bicset>0)));
for i = bicbesti
    betahat = lare(x,y,lambda(i),6);
    betaset(:,i) = betahat;
end
finalbeta = betaset(:,bicbesti);
plot(bicset)
toc
