% This Demo shows how to use functions fprcal to compute the true positive
% rate. 
% Yangguang Zang <yangguang.zang@gmail.com>
% $Revision: 1.0.0 $  $Date: 2016/05/03 $
betaori = [1;2;3;0;0;0;0;0;0];
betahat = [0.9;2.1;3.2;0.6;0.2;0.1;0;0;0];
result = fprcal(betahat,betaori);

