% This Demo shows how to use function km to estimate the weight when the 
% data is right censoring.
% Yangguang Zang <yangguang.zang@gmail.com>
% $Revision: 1.0.0 $  $Date: 2016/05/03 $
clear
n = 200;
p1 = 500;
p2 = 5;
w = randn(n,p1);
x1 = w;
x2 = randn(n,p2);
x3 = zeros(n,p1*p2);
for i = 1:n
    x3(i,:) = reshape(x2(i,:)'*x1(i,:),p1*p2,1);
end
x = [x1,x2,x3];
[~,p] = size(x);
beta01 = zeros(p1,1);
beta02 = zeros(p2,1);
beta03 = zeros(p1*p2,1);
beta01(1:10) = 0.4+0.8*rand(10,1);
beta02(1:5) = 4+0.8*rand(5,1);
beta03(1:20) = 4+0.8*rand(20,1);
beta = [beta01;beta02;beta03];
YY = exp(x*beta).*exp(randn(n,1));
C = 100*rand(n,1);
y = zeros(n,1);
Delta = zeros(n,1);
for i = 1:n
    y(i) = min(C(i),YY(i));
    if YY(i) <= C(i)
        Delta(i) = 1;
    end
end
weight = km(y,Delta);
