% NURBS-Divergence-Meshless (NDM) Method: NDM-MLS
% Title: Codes for 2D Quarter-annular Disc Heat Transfer Problem
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
mx = 8;         % Number of division in x-direction
my = 8;         % Number of division in y-direction
ng = 3;         % Number of Gpoints

% --------------------
% Prepare nodes & loops 
% --------------------

% Get NURBS parameters
ri = 0.3;
ro = 0.6;
[cpw,p,q,xi,eta] = NURBS_Qdisc(ri,ro);

% Create Nodes coordinates [x y] para-space
[x,y] = meshgrid(0:1/mx:1,0:1/my:1);
node_ab = [x(:) y(:)];
node_ab = sortrows(node_ab,2);

% Get Gauss points & associated condition
[S,C,P,W] = NDM_num(1,1,mx+1,my+1,ng);


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
Plot('Node',P(unique(S(:)),:),'bx',4)
Plot('Node',node_ab,'k.',30)

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
Plot('Node',P2(unique(S(:)),:),'bx',4)
Plot('Node',cpw(:,1:2),'r.',40)
Plot('Node',node_xy,'k.',25)


%%
% --------------------
% Calculate the local matrix
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


%% Additional plotting
figure

subplot(2,2,1) %physc-space + cont.points
axis equal tight
axis([0,0.65,0,0.65])
hold on
plot(NaN,'k.','MarkerSize',30)
plot(NaN,'b:','LineWidth',1.5)
hold off
Plot('KnotR',Obj,'r-.',3)
Plot('Node',cpw(:,1:2),'k.',40)
cpx = reshape(cpw(:,1),3,[]);
cpy = reshape(cpw(:,2),3,[]);
hold on
plot(cpx,cpy,'b:','LineWidth',1.5)
plot([],'k.','MarkerSize',10)
plot(cpx',cpy','b:','LineWidth',1.5)
hold off
legend('Control points','Control mesh')
set(gca,'FontSize',15)

subplot(2,2,2) %para-space
axis equal tight
axis([0,1.1,0,1.1])
Plot('Knot',Obj,'r-.',3)
set(gca,'FontSize',15)

subplot(2,2,3:4) %physc-space
axis equal tight off
Obj.T = T*0;
Plot('ContR',Obj)
Plot('BoundR',Obj,'k-')
set(gca,'FontSize',15)
colormap gray
clim([-30 10])
colorbar off


