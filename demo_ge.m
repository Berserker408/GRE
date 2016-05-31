% This Demo shows how to use functions lpre/lpre to select the important 
% genes, environment factors, and gene-environment interactions under the 
% AFT model. Note that in the Demo the reponse is also right censoring, so 
% we weight samples with the Kaplan-Meier weight. The best lambda is chosen
% as follows: we generate a same structure data without censoring as the test
% set, and select the tuning parameter which makes the sum of squares errors 
% of the test set smallest. We also can use the 5-fold cross validation to
% select the best tuning parameter.
% Yangguang Zang <yangguang.zang@gmail.com>
% $Revision: 1.0.0 $  $Date: 2016/05/03 $
tic,
clear
n = 200;
p1 = 50;
p2 = 5;
maxit = 20;
maxit1 = 20;
toler = 1E-4;
nLambda = 20;
rr = 6;
lambdaRatio = 1E-5;
%%%%%%%%%%%%%%%%%%%%%%%%
w = randn(n,p1);
wtest = randn(n,p1);
x1 = w;
x2 = randn(n,p2);
x3 = zeros(n,p1*p2);
x1test = wtest;
x2test = randn(n,p2);
x3test = zeros(n,p1*p2);
%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:n
    x3(i,:) = reshape(x2(i,:)'*x1(i,:),p1*p2,1);
end
for i = 1:n
    x3test(i,:) = reshape(x2(i,:)'*x1(i,:),p1*p2,1);
end
x = [x1,x2,x3];
xtest = [x1test,x2test,x3test];
[~,p] = size(x);
%%%%%%%%%%%%%%%%%%%%%%%%%
beta01 = zeros(p1,1);
beta02 = zeros(p2,1);
beta03 = zeros(p1*p2,1);
beta01(1:10) = 0.4+0.8*rand(10,1);
beta02(1:5) = 0.4+0.8*rand(5,1);
beta03(1:20) = 0.4+0.8*rand(20,1);
beta = [beta01;beta02;beta03];
%%%%%%%%%%%%%%%%%%%%%%%%%%%
YY = exp(x*beta).*exp(randn(n,1));
YYtest = exp(xtest*beta).*exp(randn(n,1));
% YY = exp(x*beta).*exp(4*rand(n,1)-2);
% epsilon = zeros(n,1);
% for i = 1:n
%     cc = rand();
%     if cc<0.9
%         epsilon(i) = randn();
%     else
%         epsilon(i) = trnd(1);
%     end
% end
% YY = exp(x*beta).*exp(epsilon);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C = 100.*rand(n,1);
y = zeros(n,1);
Delta = zeros(n,1);
for i = 1:n
    y(i) = min(C(i),YY(i));
    if YY(i) <= C(i)
        Delta(i) = 1;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yold = y;
[y,index] = sort(y);
x = x(index,:);
Delta = Delta(index,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
weight = km(y,Delta);
%%%%%%%%%%%%%%%%%%%%%%%%
r = 6;
beta0 = zeros(p,1);
for i = 1:100
    lam = i;
    betahat = lpre(x,y,lam,6);
    if sum(abs(betahat) <= 0.001)  ==  p
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
tprset = zeros(nLambda,1);
fprset = zeros(nLambda,1);
%%%%%%%%%%%%%%%%%%%%%%%%%
ebicset = zeros(nLambda,1);
for i = 1:nLambda
    i
	beta0 = lpre(x,y,lambda(i),6);
    betaset(:,i) = beta0;
    tprset(i,1) = tprcal(beta0,beta);
    fprset(i,1) = fprcal(beta0,beta);
    ebicset(i,1) = 2*log(sum((log(YYtest)-xtest*beta0).^2));
end
auc = auccalc(fprset,tprset);
ebicbesti = find(ebicset == min(ebicset(ebicset>0)));
finalbeta = betaset(:,ebicbesti);
plot(ebicset)
toc

