% Title: Calculate the local matrix and force vector for Cellbeam unit
% Created by: Malaysian Society for Numerical Methods (MSNM)
% Copyright © 2024 . All Rights Reserved


function [K,F,node_xy] = Calc_kf_UnitCB(kx,ky,qt,ng,mx,my,cpw,p,q,xi,eta)

% --------------------
% Calculate the local matrix
% --------------------

% Create Nodes coordinates [x y] para-space
[x,y] = meshgrid(0:1/mx:1,0:1/my:1);
node_ab = [x(:) y(:)];
node_ab = sortrows(node_ab,2);

% Get Gauss points & associated condition
[S,C,P,W] = NDM_num(1,1,mx+1,my+1,ng);

% Get shape function to calculate dNx, dNf
h = max([1/mx 1/my]);
rc = 5*h;
[~,dNdx,dNdy] = MLS_SF(P,node_ab,rc);

%calculate for G1 & G2 matrix
[~,dRds,dRdn] = NURBS_SF(P,cpw,p,q,xi,eta);
dxds = dRds*cpw(:,1);
dxdn = dRdn*cpw(:,1);
dyds = dRds*cpw(:,2);
dydn = dRdn*cpw(:,2); 
detJ = abs(dxds.*dydn - dxdn.*dyds);

G11 = (dxdn.^2 + dydn.^2) ./ detJ;
G12 = -(dxds.*dxdn + dyds.*dydn) ./ detJ;
G22 = (dxds.^2 + dyds.^2) ./ detJ;

G1 = [-G12 G11 G12 -G11];
G2 = [-G22 G12 G22 -G12];

%evaluated dNx & dNy at P
dNx = kx*dNdx;
dNy = ky*dNdy;

K = zeros(size(S,1),size(dNdx,2));
F = zeros(size(S,1),1);
for i = 1:size(S,1) %loop for sub-area
    s = S(i,:);     %related sides for loop i

    % Generate G1 & G2 matrix for a loop 
    g1 = []; g2 = [];
    for j = 1:4
        ind = ng*(j-1)+1;
        sj = s(ind:ind+ng-1);
        g1 = cat(1,g1,G1(sj,j));
        g2 = cat(1,g2,G2(sj,j));
    end

    % Calculate K matrix
    for j = 1:size(dNdx,2) %loop over SF
        
        % combine all derivaties for each side
        dNxi = sum(reshape(W(s,1).*dNx(s,j).*g1,ng,[]),1)';
        dNyi = sum(reshape(W(s,1).*dNy(s,j).*g2,ng,[]),1)';
        
        % remove boundaries term
        dNxi = dNxi.*C(i,:)';
        dNyi = dNyi.*C(i,:)';
        
        % sum up & put in k matrix
        dN = sum(dNxi) + sum(dNyi);
        K(i,j) = dN;
    end
    
end

% ----------------------
% Impose boundary conditions (Neumann)
% ----------------------

for i = 1:size(S,1) %loop for sub-area
    s = S(i,:);     %related sides for loop i
    bt = P(s,2)==max(P(:,2));   %index for top boundary points
    dRt = dRds(s,:);    %derivative SF for top boundary points
    detJi = sqrt((dRt*cpw(:,1)).^2+(dRt*cpw(:,2)).^2);
    F(i,1) = F(i,1) -sum(W(s,1).*bt.*detJi*qt);
end

R1 = NURBS_SF(node_ab,cpw,p,q,xi,eta);
node_xy = R1*cpw(:, 1:2);


end