% NURBS-Divergence-Meshless (NDM) Method: NDM-MLS
% Title: Codes for 2D Fan Blade Heat Transfer Problem - Case 2&3
% Created by: Malaysian Society for Numerical Methods (MSNM)
% Copyright © 2024 . All Rights Reserved

clc
clear
close all

% ----------------------
% Input: 
% ----------------------

% The parameters
kx = 1;         % Thermal conductivity in x direction [W/(m.degC)]
ky = 1;         % Thermal conductivity in y direction [W/(m.degC)]
Q  = 0;         % Rate of heat generation [W/m^3]
ql = 60;        % Heat flux for left boundary [W/m^2]
Tr = 30;        % T along right boundary [degC]
Tt = 30;        % T along top boundary [degC]
Tb = 30;        % T along bottom boundary [degC]
mx = 12;        % Number of division in x-direction
my = 12;        % Number of division in y-direction
ng = 1;         % Number of Gpoints

% --------------------
% Prepare nodes & looping 
% --------------------

% Get NURBS parameters
[cpw,p,q,xi,eta] = NURBS_FanBlade();

% Create Nodes coordinates [x y] para-space
[x,y] = meshgrid(0:1/mx:1,0:1/my:1);
node_ab = [x(:) y(:)];
node_ab = sortrows(node_ab,2);

% Get Gauss points & associated condition
[S,C,P,W] = NDM_num(1,1,mx+1,my+1,ng);

% Case 2 - varies location of nodes
n1 = [];
for ON = [ 1 ]
    n1 = [29 30 31 37 42 48 49 74 75 76 94 95 113 127 134 135 141];
    dx = [-7  1  5  5 -5 -6 7  -8 -8 -7 -6 -7   7   8  -8   4  -5]/my/20;
    dy = [-8 -11 -9 8 5 -11 8 -12 18 -13 -7 13 -4  -8   7   8   9]/my/30;
    for i = 1:length(n1)
        j = n1(i);
        node_ab(j,:) = node_ab(j,:) + [dx(i) dy(i)];
    end
end

% Case 3 - varies size of loops
n2 = [];
for ON = [ 1 ]
    Smax = max(S(:));
    n2 = [15 18 35 24 70 71 100 102 115 68 146 148 139 149 86 86 45];
    dxy = abs(dx)+abs(dy)/1.7;
    for i = 1:length(n2)
        pm = P(S(n2(i),:),:) + [0 -1; 1 0; 0 1; -1 0]*dxy(i);
        P = cat(1,P,pm);
        wm = [pm(2,1)-pm(4,1)
            pm(3,2)-pm(1,2)
            pm(2,1)-pm(4,1)
            pm(3,2)-pm(1,2)];
        W = cat(1,W,wm);
        S(n2(i),:) = (1:4*ng)+Smax;
        Smax = max(S(:));
    end
    newnode_ab = node_ab(n2,:);
end


% Plotting
Obj.xi = xi;
Obj.eta = eta;
Obj.S = S;
Obj.P = P;
Obj.cpw = cpw;
Obj.p = p;
Obj.q = q;

% parametric space
figure
subplot(1,2,1)
axis equal tight
axis([0,1,0,1])
Plot('Knot',Obj,'r-.',3)
Obj.S = S;
Plot('Elem',Obj,'b--')
Obj.S = S(n2,:);
Plot('Elem',Obj,'m--')
Plot('Node',P(unique(S(:)),:),'bx',4)
Plot('Node',node_ab,'k.',30)
Plot('Node',node_ab(n1,:),'m.',30)

% physical space
R1 = NURBS_SF(node_ab,cpw,p,q,xi,eta);
node_xy = R1*cpw(:, 1:2);
R2 = NURBS_SF(P,cpw,p,q,xi,eta);
P2 = R2*cpw(:, 1:2);

subplot(1,2,2)
axis equal tight
Plot('KnotR',Obj,'r-.',3)
Obj.S = S;
Plot('ElemR',Obj,'b--')
Obj.S = S(n2,:);
Plot('ElemR',Obj,'m--')
Plot('Node',P2(unique(S(:)),:),'bx',4)
Plot('Node',cpw(:,1:2),'r.',40)
Plot('Node',node_xy,'k.',25)
Plot('Node',node_xy(n1,:),'m.',30)



%%
% --------------------
% Calculate the local matrix and force vector
% --------------------

% Get shape function to calculate dNx, dNf
h = max([1/mx 1/my]);
ac = 1;
dc = h;
rc = 4*h;
[N,dNdx,dNdy] = MLS_SF(P,node_ab,rc);

%calculate for G1 & G2 matrix
[Ri,dRds,dRdn] = NURBS_SF(P,cpw,p,q,xi,eta);
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

Line = 0;
for i = 1:size(S,1) %loop for sub-area
    s = S(i,:);     %related sides for loop i
    disp([i s])
    bl = P(s,1)==0;   %index for left boundary points
    dRl = dRdn(s,:);  %derivative SF for left boundary points
    detJi = sqrt((dRl*cpw(:,1)).^2+(dRl*cpw(:,2)).^2);
    F(i,1) = F(i,1) -sum(W(s,1).*bl.*detJi*ql);
end

% ----------------------
% Impose boundary conditions (Drichlet)
% ----------------------

% Nodes boundary conditions [Condition(1-boundary node, 0-node to be solved)  Value]
Rxy = NURBS_SF(node_ab,cpw,p,q,xi,eta);
node_xy = Rxy*cpw(:,1:2);
node_bc = zeros(length(node_ab),2);
node_bc(node_ab(:,1)==1,1) = 1;
node_bc(node_ab(:,1)==1,2) = Tr;
node_bc(node_ab(:,2)==0,1) = 1;
node_bc(node_ab(:,2)==0,2) = Tb;
node_bc(node_ab(:,2)==1,1) = 1;
node_bc(node_ab(:,2)==1,2) = Tt;

for i = 1:size(node_ab,1)
    % Check if node condition is boundary
    if node_bc(i,1) > 0       
        K(i,:) = 0;
        K(i,i) = 1;
        F(i)   = node_bc(i,2);        
    end
end

% ----------------------
% Solve simultaneous equations 
% ----------------------

T = K\F;
T0 = T(1:(mx+1)*(my+1));
T1 = reshape(T0,mx+1,[])'
fprintf('max T = %.3f degC\n',max(T))

% --------------------
% Plot the result
% --------------------

% close all
figure('color', 'w', 'Name', 'Plot of Temperature Distribution');
axis equal;

Obj.mx = mx;
Obj.my = my;
Obj.node_ab = node_ab;
Obj.T = T;
Plot('ContR',Obj)
Plot('KnotR',Obj,'r-',3)

Obj.x = 0:1/mx:1; Obj.y = 0:1/my:1; 
Plot('Node',node_xy,'k.',15)
xlabel('x [m]')
ylabel('y [m]')

