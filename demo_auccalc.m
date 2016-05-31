% This Demo shows how to use functions auccalc to compute the area under
% the ROC (AUC).
% Yangguang Zang <yangguang.zang@gmail.com>
% $Revision: 1.0.0 $  $Date: 2016/05/03 $
fprset=[0; 0; 0.0111;0.0444;0.1037;
    0.2407;0.3778;0.4519;0.5111;0.5667;
    0.6000;0.6185; 0.6481;0.6519;0.6889;
    0.7222;0.7889;0.8296;0.9481;0.9852];
tprset=[0;0.0286;0.3429;0.6571;0.8000;
    0.8857;0.8857;0.8571;0.9143;0.9143;
    0.9143;0.9143;0.9143;0.9143;0.9143;
    0.9143;0.9714;1.0000;1.0000;1.0000];
auc = auccalc(fprset,tprset);

