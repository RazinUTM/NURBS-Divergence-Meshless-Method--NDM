% Title: Quarter Disc benchmarking problem
% Attempt: to make a module - simplify the main code
% Created on: 13/08/2021 by Razin
% Last updated: 13/08/2021 by Razin
% Double checked on: 13/08/2021 by Razin


function [cpw,p,q,xi,eta] = NURBS_Qdisc(ri,ro)


% Input NURBS
cpw = [ri                          0              1;
       ri+(ro-ri)/2                0              1;
       ro                          0              1;
       ri                          ri*(sqrt(2)-1) 0.5*(1+1/sqrt(2));
       ri+(ro-ri)/2  (ri+(ro-ri)/2)*(sqrt(2)-1)   0.5*(1+1/sqrt(2));
       ro                          ro*(sqrt(2)-1) 0.5*(1+1/sqrt(2));
       ri*(sqrt(2)-1)              ri             0.5*(1+1/sqrt(2));
       (ri+(ro-ri)/2)*(sqrt(2)-1)  ri+(ro-ri)/2   0.5*(1+1/sqrt(2));
       ro*(sqrt(2)-1)              ro             0.5*(1+1/sqrt(2));
       0                           ri             1;
       0                           ri+(ro-ri)/2   1;
       0                           ro             1];   
p = 2;                       %poly order in x-dir
q = 2;                       %poly order in y-dir
xi = [0 0 0 1 1 1];          %knot in x-dir (natural space)
eta = [0 0 0 0.5 1 1 1];     %knot in y-dir (natural space)

end