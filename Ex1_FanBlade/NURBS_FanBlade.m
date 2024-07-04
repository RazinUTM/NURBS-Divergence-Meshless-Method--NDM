% Title: Fan Blade benchmarking problem - Luo et.al. (2023)
% Created by: Malaysian Society for Numerical Methods (MSNM)
% Copyright © 2024 . All Rights Reserved


function [cpw,p,q,xi,eta] = NURBS_FanBlade()


% Input NURBS
cpw = [0.0515    0.0037    1.0000
       0.0937    0.0014    0.9821
       0.1320   -0.0163    1.0000
       0.1868   -0.0417    0.4984
       0.1799    0.0271    1.0000
       0.0470    0.0651    0.6421
       0.0883    0.1015    0.6796
       0.1296    0.1378    0.7170
       0.1466    0.1300    0.8157
       0.1636    0.1222    0.9143
      -0.0127    0.0500    1.0000
      -0.0207    0.0956    0.9882
      -0.0085    0.1478    1.0000
       0.0305    0.2048    0.6698
       0.0939    0.1661    1.0000]; 
  
p = 2;                           %poly order in x-dir
q = 2;                           %poly order in y-dir
xi = [0 0 0 0.5 0.5 1 1 1];      %knot in x-dir (natural space)
eta = [0 0 0 1 1 1];             %knot in y-dir (natural space)

end
