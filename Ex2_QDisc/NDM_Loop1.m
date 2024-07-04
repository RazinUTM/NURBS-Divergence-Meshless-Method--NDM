% Title: Create One Loop (FVM)
% Attempt: to make a loop with any Gpoints
% Created on: 16/12/2023 by Razin
% Last updated: 16/12/2023 by Razin


function [P,W] = NDM_Loop1(Nod,L,ng)

% function [P,W] = NDM_Loop1()
% Nod = [0 0];
% L = [1 2];
% ng = 3;


% Input
x0 = Nod(1);
y0 = Nod(2);
lx = L(1);
ly = L(2);

% Integration points to evaluate
[xg,wg] = Gauss_Points_1D(ng,1);

o1 = ones(ng,1);
z0 = zeros(ng,1);
%generate the location of points
xgb = [ xg*lx  z0   ];
xgr = [ o1*lx  xg*ly ];
xgt = [ xg(end:-1:1)*lx  o1*ly ];
xgl = [ z0     xg(end:-1:1)*ly ];

P0 = [xgb; xgr; xgt; xgl];
P = [P0(:,1)-lx/2+x0 , P0(:,2)-ly/2+y0];
% AddPlot2(P,'rx',12)

%adjust weight Gpoints
wgb = wg*lx;
wgr = wg*ly;
wgt = wg(end:-1:1)*lx;
wgl = wg(end:-1:1)*ly;
W = [wgb; wgr; wgt; wgl];

end

function [xg,wg] = Gauss_Points_1D(ng,L)
% This function
% calculates Gpoints [xg] & weight [wg] for 1D
% for a given no. of Gpoints ng x ng (ng)
% and the line length (L)

c = (1:ng-1)./sqrt(4*(1:ng-1).^2-1);
J = diag(c,-1) + diag(c,1); % Diagonal values is 0
[V,D] = eig(J);
[X,ix] = sort(diag(D)); % Sort diagonal eigenvalues in ascending order for x

V(:,1:ng) = V(:,ix);
w = 2*V(1,:)'.^2;

xg = X*(L/2) + (L/2); % change of interval
wg = w*L/2;

end