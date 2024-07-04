% NURBS-Divergence-Meshless (NDM) Method: NDM-MLS
% Title: Codes for 2D Cellular Web Problem
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
qt = 60;        % Heat flux for top side [W/m^2]
Tb = 30;        % T along bottom side [degC]
Tl = 30;        % T along left side [degC]
Tr = 30;        % T along right side [degC]
mx = 20;        % Number of division in x-direction - per unit
my = 20;        % Number of division in y-direction
ng = 3;         % Number of Gpoints

% Cellbeam dimension
d = 0.57*2;     % Depth [m]
r = 0.4;        % Radius [m]
sL = 1.06;      % Width - Left CB Unit [m]
sM = 1.47;      % Width - Mid CB Unit [m]
sR = 1.06;      % Width - Right CB Unit [m]

% --------------------
% Prepare nodes & loops,
% Calculate the local stiffness matrix and force vector
% --------------------

% UnitCB - Left
[cpw,p,q,xi,eta] = NURBS_UnitCB_Left(sL,d,r);
[KL,FL,node_xyL] = Calc_kf_UnitCB(kx,ky,qt,ng,mx,my,cpw,p,q,xi,eta);

% UnitCB - Mid
[cpw,p,q,xi,eta] = NURBS_UnitCB_Mid(sM,d,r);
[KM,FM,node_xyM] = Calc_kf_UnitCB(kx,ky,qt,ng,mx,my,cpw,p,q,xi,eta);

% UnitCB - Right
[cpw,p,q,xi,eta] = NURBS_UnitCB_Right(sR,d,r);
[KR,FR,node_xyR] = Calc_kf_UnitCB(kx,ky,qt,ng,mx,my,cpw,p,q,xi,eta);

% --------------------
% Obtain the global stiffness matrix and force vector
% --------------------

ShiftX = [0 sL sL+sM sL+2*sM sL+3*sM]; % adjust UnitCB location
node_xy_G = [];
Klocal = [];
Flocal = [];
for k = 1:length(ShiftX)
    
    if k == 1
        node_xy = node_xyL;
        Klocal(k).Klocal = KL;
        Flocal(k).Flocal = FL;
        
    elseif k == length(ShiftX)
        node_xy = node_xyR;
        Klocal(k).Klocal = KR;
        Flocal(k).Flocal = FR;
        
    else
        node_xy = node_xyM;
        Klocal(k).Klocal = KM;
        Flocal(k).Flocal = FM;
    end
    
    % Shift coordinate x accordingly
    node_xy(:,1) = node_xy(:,1) + ShiftX(k)*ones(size(node_xy(:,1),1),1);
    node_xy_G = cat(1,node_xy_G,node_xy);

end

node_xy_G0 = node_xy_G;
node_xy_G = round(node_xy_G,5);
[node_xy_G,~,iori] = unique(node_xy_G,'rows');
TN = size(node_xy_G,1);     % Total nodes

% ----------------------
% Assemble into Global
% ----------------------
num = 1:size(node_xy_G,1);  % list 1... Total nodes
num = num(iori);            % connectivity for nodes 
numP = num;

FGlobal = zeros(TN,1);
KGlobal = zeros(TN,TN);
for i=1:length(ShiftX)
    
    % Extract nodes connectivity for each patch
    n1 = size(Klocal(i).Klocal,1);
    uv = num(1:n1); % nodes connectivity
    
    % Assemble into Global
    FGlobal(uv,:)   = FGlobal(uv,:) + Flocal(i).Flocal;
    KGlobal(uv,uv)  = KGlobal(uv,uv)+ Klocal(i).Klocal;
    num(1:n1) = [];
    
end

% ----------------------
% Impose boundary conditions (Drichlet)
% ----------------------
TNodeList = [];
for iNode = 1:TN
    
    x = node_xy_G(iNode,1);
    y = node_xy_G(iNode,2);

    if (y==0 || x==0 || x>=max(node_xy(:,1))-1e-3) %&& y<=max(node_xy(:,2))-1e-3
        
        KGlobal(iNode,:)	= 0;
        KGlobal(iNode,iNode) = 1;
        FGlobal(iNode)      = Tb;
    end
    
end

% ----------------------
% Solve simultaneous equations 
% ----------------------
T = KGlobal\FGlobal;

% ----------------------
% Display result at selected points
% ----------------------
node_xy_G1 = round(node_xy_G,3);
fprintf('   x [m]    y [m]     T ["C]\n\n')
for iNode = 1:TN

    x = node_xy_G1(iNode,1);
    y = node_xy_G1(iNode,2);

    if (y==0.97 && x==2.53) || (y==1.14 && x==2.53) || ...
            (y==0.57 && x==2.93) || (y==1.14 && x==3.265)

        Ti = T(iNode,1);
        fprintf('%8.3f %8.3f %10.3f\n',x,y,Ti)
    end

end

%% Plotting 01 - Nodes

figure('color', 'w', 'Name', 'Plot Nodes');
hold on
for i = 1:length(ShiftX)
    
    if i == 1
        node_xyP = node_xyL;
        [cpw,p,q,xi,eta] = NURBS_UnitCB_Left(sL,d,r);      
    elseif i == length(ShiftX)
        node_xyP = node_xyR;
        [cpw,p,q,xi,eta] = NURBS_UnitCB_Right(sR,d,r);
    else
        node_xyP = node_xyM;
        [cpw,p,q,xi,eta] = NURBS_UnitCB_Mid(sM,d,r);
    end
    [S,C,P,W] = NDM_num(1,1,mx+1,my+1,ng);

    % Shift coordinate x accordingly
    ShiftXi = ShiftX(i);
    node_xyP(:,1) = node_xyP(:,1) + ShiftXi;
    node_xy = node_xyP;

    % required variables
    Obj.xi = xi;
    Obj.eta = eta;
    Obj.cpw = cpw;
    Obj.cpw(:,1) = Obj.cpw(:,1) + ShiftXi;
    Obj.p = p;
    Obj.q = q;
    Obj.S = S;
    Obj.P = P;

    axis([0,1.5,0,1.2])
    Plot('ElemR',Obj,'b--')
    Plot('KnotR',Obj,'r',2)

end
Plot('Node',node_xy_G,'k.',10)
hold off
set(gcf, 'Position',  [100, 100, 1100, 400])
set(gca,'FontSize',12)
axis tight equal


%% Figure 02 - Contour plot

figure('color', 'w', 'Name', 'Plot Contour');
numPi = numP;
hold on
for k = 1:length(ShiftX)
    
    if k == 1
        node_xyP = node_xyL;
        [cpw,p,q,xi,eta] = NURBS_UnitCB_Left(sL,d,r);        
    elseif k == length(ShiftX)
        node_xyP = node_xyR;
        [cpw,p,q,xi,eta] = NURBS_UnitCB_Right(sR,d,r);
    else
        node_xyP = node_xyM;
        [cpw,p,q,xi,eta] = NURBS_UnitCB_Mid(sM,d,r);
    end
    [S,C,P,W] = NDM_num(1,1,mx+1,my+1,ng);
    
    % Shift coordinate x accordingly
    node_xyP(:,1) = node_xyP(:,1) + ShiftX(k)*ones(size(node_xyP(:,1),1),1);

    % Extract T
    n1 = size(Klocal(k).Klocal,1);
    uv = numPi(1:n1); % nodes connectivity
    TP = T(uv,1);
    numPi(1:n1) = [];

    % Nodes coordinates [a b]
    [a,b] = meshgrid(0:1/mx:1,0:1/my:1);
    node_ab = [a(:) b(:)];
    node_ab = sortrows(node_ab,2);

    % Nodes coordinates [x y]
    plot_mx = 40;
    plot_my = plot_mx;
    [x,y] = meshgrid(0:1/plot_mx:1,0:1/plot_my:1);
    plot_ab = [x(:) y(:)];
    plot_ab = sortrows(plot_ab,2);
    R3 = NURBS_SF(plot_ab,cpw,p,q,xi,eta);
    plot_xy = R3*cpw(:, 1:2);
    plot_xy(:,1) = plot_xy(:,1) + ShiftX(k)*ones(size(plot_xy(:,1),1),1);

    % Interpolation
    N = MLS_SF(plot_ab,node_ab,0.6);
    plot_T = N*TP;

    % Plot 2D contour of the result
    elem_node = [];
    for i = 1:plot_mx
        for j = 1:plot_mx
            elem_i = [[j j+1]  [j+1 j]+(plot_mx+1)] + (plot_mx+1)*(i-1);
            elem_node = cat(1,elem_node,elem_i);
        end
    end

    patch('Faces', elem_node, ...
        'Vertices', plot_xy, ...
        'FaceVertexCData', plot_T, ...
        'EdgeColor','none', ...
        'FaceColor','interp');
    colormap('jet')
    colorbar

end
hold off
set(gcf, 'Position',  [100, 100, 1100, 400])
set(gca,'FontSize',12)
axis tight equal
