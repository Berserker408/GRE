function weight=km(y,Delta)
% Function km creates the Kaplan-Meier weight for right censoring data.
% Input:
% 	y: n dimensional right censoring observation.
%	Delta:  n dimensional event indicator, which has the same length with y. The
%           element of Delta equals 1 if the event was observed, or equals 0 if
%           the responce was censored
% Ouput:
%   weight: n dimensional Kaplan-Meier vector.
% Yangguang Zang <yangguang.zang@gmail.com>
% $Revision: 1.0.0 $  $Date: 2016/05/03 $
n=size(y,1);
[~,index]=sort(y);
Delta=Delta(index,:);
weight=zeros(n,1);
for i=1:n
    if i==1
        weight(i)=Delta(i)/n;
    else
        temp=zeros(i-1,1);
        for j=1:i-1
            temp(j)=((n-j)/(n-j+1))^Delta(j);
        end
        weight(i)=Delta(i)/(n-i+1)*prod(temp);
    end
end