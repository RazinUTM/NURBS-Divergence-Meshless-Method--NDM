% Title: Plot Nodes & Elements - centralised all plotting functions
% Created by: Malaysian Society for Numerical Methods (MSNM)
% Copyright © 2024 . All Rights Reserved


function Plot(Opt,Obj,Prop,Scale)

hold on
switch Opt
    case 'Knot'
        Knot(Obj,Prop,Scale)   % plot Knot spans
    case 'Node'
        Node(Obj,Prop,Scale)   % plot nodes/points
    case 'Text'
        Text(Obj,Prop)         % plot nodal number
    case 'Elem'
        Elem(Obj,Prop)         % plot elements
    case 'Elem2'
        Elem2(Obj,Prop)        % plot elements (with numbers)

    case 'KnotR'
        KnotR(Obj,Prop,Scale)  % plot Knot spans Nurbs
    case 'BoundR'
        BoundR(Obj,Prop)       % plot boundary Nurbs
    case 'MeshR'
        MeshR(Obj,Prop)        % plot Knot spans Nurbs
    case 'NodeR'
        NodeR(Obj,Prop,Scale)  % plot nodes/points Nurbs
    case 'ElemR'
        ElemR(Obj,Prop)        % plot elements Nurbs
    case 'ElemR2'
        ElemR2(Obj,Prop)       % plot elements Nurbs (with numbers)
    case 'ContR'
        ContR(Obj)             % plot contour of elements Nurbs
end
hold off

end

function Knot(Obj,Prop,Scale)
% plot Knot spans
xi = Obj.xi;
eta = Obj.eta;
Xi = unique(xi);
Eta = unique(eta);

Lw = Scale; %linewidth
Ls =@(s1,s2) linspace(s1,s2,10);
for i = 1:length(Xi)-1
    for j = 1:length(Eta)

        X = Ls(Xi(i),Xi(i+1))';
        Y = X*0+Eta(j);
        XY = [X Y];

        plot(XY(:,1),XY(:,2),Prop,'LineWidth',Lw)

    end
end
for i = 1:length(Eta)-1
    for j = 1:length(Xi)

        Y = Ls(Eta(i),Eta(i+1))';
        X = Y*0+Xi(j);
        XY = [X Y];

        plot(XY(:,1),XY(:,2),Prop,'LineWidth',Lw)

    end
end
end

function Node(Obj,Prop,Scale)
% plot Nodes / Points
plot(Obj(:,1),Obj(:,2),Prop,'markerSize',Scale)
end

function Text(Obj,Prop)
% plot nodal number
for i = 1:size(Obj,1) %all integration points
    x = Obj(i,1);
    y = Obj(i,2);
    text(x,y,num2str(i),'Color',Prop,'fontsize',10, ...
        VerticalAlignment="top", HorizontalAlignment="left")
    % pause(0.1)
end
end

function Elem(Obj,Prop)
% plot elements
P = Obj.P;  S = Obj.S;
for i = 1:size(S,1) %elements (sub-area)
    x1 = min(P(S(i,:),1));
    x2 = max(P(S(i,:),1));
    y1 = min(P(S(i,:),2));
    y2 = max(P(S(i,:),2));
    X = [x1 x2 x2 x1 x1];
    Y = [y1 y1 y2 y2 y1];
    plot(X,Y,Prop,'LineWidth',0.1)
end
end

function Elem2(Obj,Prop)
% plot elements Nurbs
P = Obj.P;  S = Obj.S;
for i = 1:size(S,1) %elements (sub-area)
    x1 = min(P(S(i,:),1));
    x2 = max(P(S(i,:),1));
    y1 = min(P(S(i,:),2));
    y2 = max(P(S(i,:),2));
    X = [x1 x2 x2 x1 x1];
    Y = [y1 y1 y2 y2 y1];
    plot(X,Y,Prop,'LineWidth',0.1)

    x = (x1+x2)/2;
    y = (y1+y2)/2;
    text(x,y,num2str(i),'Color',[0.7, 0.8 ,1],'fontsize',10, ...
        VerticalAlignment="middle", HorizontalAlignment="center")
    % pause(0.1)
end
end

function KnotR(Obj,Prop,Scale)
% plot Knot spans
cpw = Obj.cpw;
p = Obj.p;
q = Obj.q;
xi = Obj.xi;
eta = Obj.eta;
Xi = unique(xi);
Eta = unique(eta);

Lw = Scale; %linewidth
Ls =@(s1,s2) linspace(s1,s2,25);
for i = 1:length(Xi)-1
    for j = 1:length(Eta)

        X = Ls(Xi(i),Xi(i+1))';
        Y = X*0+Eta(j);

        Rx = NURBS_SF([X Y],cpw,p,q,xi,eta);
        XY = Rx*cpw(:, 1:2);

        plot(XY(:,1),XY(:,2),Prop,'LineWidth',Lw)

    end
end
for i = 1:length(Eta)-1
    for j = 1:length(Xi)

        Y = Ls(Eta(i),Eta(i+1))';
        X = Y*0+Xi(j);

        Rx = NURBS_SF([X Y],cpw,p,q,xi,eta);
        XY = Rx*cpw(:, 1:2);

        plot(XY(:,1),XY(:,2),Prop,'LineWidth',Lw)

    end
end
end

function BoundR(Obj,Prop)
% plot Knot spans
cpw = Obj.cpw;
p = Obj.p;
q = Obj.q;
xi = Obj.xi;
eta = Obj.eta;
Xi = unique(xi); Xi = Xi([1 end]);
Eta = unique(eta); Eta = Eta([1 end]);

Lw = 1;
Ls =@(s1,s2) linspace(s1,s2,100);
for i = 1:length(Xi)-1
    for j = 1:length(Eta)

        X = Ls(Xi(i),Xi(i+1))';
        Y = X*0+Eta(j);

        Rx = NURBS_SF([X Y],cpw,p,q,xi,eta);
        XY = Rx*cpw(:, 1:2);

        plot(XY(:,1),XY(:,2),Prop,'LineWidth',Lw)

    end
end
for i = 1:length(Eta)-1
    for j = 1:length(Xi)

        Y = Ls(Eta(i),Eta(i+1))';
        X = Y*0+Xi(j);

        Rx = NURBS_SF([X Y],cpw,p,q,xi,eta);
        XY = Rx*cpw(:, 1:2);

        plot(XY(:,1),XY(:,2),Prop,'LineWidth',Lw)

    end
end
end

function MeshR(Obj,Prop)
% plot Knot spans
cpw = Obj.cpw;
p = Obj.p;
q = Obj.q;
xi = Obj.xi;
eta = Obj.eta;
Xi = unique(Obj.x);
Eta = unique(Obj.y);

Ls =@(s1,s2) linspace(s1,s2,10);
for i = 1:length(Xi)-1
    for j = 1:length(Eta)

        X = Ls(Xi(i),Xi(i+1))';
        Y = X*0+Eta(j);

        Rx = NURBS_SF([X Y],cpw,p,q,xi,eta);
        XY = Rx*cpw(:, 1:2);

        plot(XY(:,1),XY(:,2),Prop,'LineWidth',0.1)

    end
end
for i = 1:length(Eta)-1
    for j = 1:length(Xi)

        Y = Ls(Eta(i),Eta(i+1))';
        X = Y*0+Xi(j);

        Rx = NURBS_SF([X Y],cpw,p,q,xi,eta);
        XY = Rx*cpw(:, 1:2);

        plot(XY(:,1),XY(:,2),Prop,'LineWidth',0.1)

    end
end
end

function NodeR(Obj,Prop,Scale)
% plot Nodes / Points
plot(Obj(:,1),Obj(:,2),Prop,'markerSize',Scale)
end

function ElemR(Obj,Prop)
% plot elements Nurbs
cpw = Obj.cpw;
p = Obj.p;
q = Obj.q;
xi = Obj.xi;
eta = Obj.eta;

Ls =@(s1,s2) linspace(s1,s2,10);

P = Obj.P;  S = Obj.S;
for i = 1:size(S,1) %elements (sub-area)

    x1 = min(P(S(i,:),1));
    x2 = max(P(S(i,:),1));
    y1 = min(P(S(i,:),2));
    y2 = max(P(S(i,:),2));
    lx1 = Ls(x1,x2);  ly1 = Ls(y1,y1);
    lx2 = Ls(x2,x2);  ly2 = Ls(y1,y2);
    lx3 = Ls(x2,x1);  ly3 = Ls(y2,y2);
    lx4 = Ls(x1,x1);  ly4 = Ls(y2,y1);
    X = [lx1 lx2 lx3 lx4]';
    Y = [ly1 ly2 ly3 ly4]';

    Rx = NURBS_SF([X Y],cpw,p,q,xi,eta);
    XY = Rx*cpw(:,1:2);

    plot(XY(:,1),XY(:,2),Prop,'LineWidth',0.1)
    % pause(0.1)
end
end

function ElemR2(Obj,Prop)
% plot elements Nurbs (with numbers)
cpw = Obj.cpw;
p = Obj.p;
q = Obj.q;
xi = Obj.xi;
eta = Obj.eta;

Ls =@(s1,s2) linspace(s1,s2,10);

P = Obj.P;  S = Obj.S;
for i = 1:size(S,1) %elements (sub-area)

    x1 = min(P(S(i,:),1));
    x2 = max(P(S(i,:),1));
    y1 = min(P(S(i,:),2));
    y2 = max(P(S(i,:),2));
    lx1 = Ls(x1,x2);  ly1 = Ls(y1,y1);
    lx2 = Ls(x2,x2);  ly2 = Ls(y1,y2);
    lx3 = Ls(x2,x1);  ly3 = Ls(y2,y2);
    lx4 = Ls(x1,x1);  ly4 = Ls(y2,y1);
    X = [lx1 lx2 lx3 lx4]';
    Y = [ly1 ly2 ly3 ly4]';
    
    Rx = NURBS_SF([X Y],cpw,p,q,xi,eta);
    XY = Rx*cpw(:,1:2);

    plot(XY(:,1),XY(:,2),Prop)

    xx = (min(XY(:,1))+max(XY(:,1)))/2;
    yy = (min(XY(:,2))+max(XY(:,2)))/2;
    text(xx,yy,num2str(i), ...
        'Color',[0.7, 0.8 ,1],'fontsize',10, ...
        VerticalAlignment="middle", HorizontalAlignment="center")
    % pause(0.1)
end
end

function ContR(Obj)
% plot contour of elements Nurbs

% get variables
cpw = Obj.cpw;
p = Obj.p;
q = Obj.q;
xi = Obj.xi;
eta = Obj.eta;
node_ab = Obj.node_ab;
T = Obj.T;

% Nodes coordinates [x y]
plot_mx = 50;
plot_my = plot_mx;
[x,y] = meshgrid(0:1/plot_mx:1,0:1/plot_my:1);
plot_ab = [x(:) y(:)];
plot_ab = sortrows(plot_ab,2);
R3 = NURBS_SF(plot_ab,cpw,p,q,xi,eta);
plot_xy = R3*cpw(:, 1:2);

% Interpolation
N = MLS_SF(plot_ab,node_ab,1);
plot_T = N*T;

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