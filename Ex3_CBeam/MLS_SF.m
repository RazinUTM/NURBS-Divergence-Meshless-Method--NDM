% Title: 2D SF Code considering MLS / Meshfree-EFG
% Created by: Malaysian Society for Numerical Methods (MSNM)
% Copyright © 2024 . All Rights Reserved

function [Nx,dNx,dNy] = MLS_SF(points,nodes,r)

%shape parameters in EFG
c = r/4;

%generate matrix P
x = nodes(:,1);
y = nodes(:,2);
p0 = zeros(size(x));
p1 = ones(size(x));
Pm = [p1 x y x.^2 x.*y y.^2]';
nP = size(Pm,1);

%calculate the shape function & derivatives
Nx = zeros(length(points),length(nodes));
dNx = zeros(length(points),length(nodes));
dNy = zeros(length(points),length(nodes));
for i = 1:length(points) %loop over gauss points
  
    %weight function, w
    d = sqrt((points(i,1)-nodes(:,1)).^2 + ...
        (points(i,2)-nodes(:,2)).^2)'; %distance
    wi = zeros(size(d));
    wi(d<=r) = (exp(-(d(d<=r)/c).^2)-...
        exp(-(r/c)^2))/(1-exp(-(r/c)^2));
        
    %derivatives wx & wy
    k = (d>0 & d<r);
    dd1 = (points(i,1)-nodes(k,1))';
    dd2 = (points(i,2)-nodes(k,2))';
    dw1 = zeros(size(d));
    dw2 = zeros(size(d));
    dw = @(i,dd) -(2*dd.*exp(-d(k).^2/c^2))/(c^2*(1-exp(-(r/c)^2)));
    dw1(k) = dw(k,dd1);
    dw2(k) = dw(k,dd2);

    %calculate matrix A & B (vectorised)
    j = wi~=0;
    P = Pm(:,j);
    Wi = (ones(nP,1)*wi(j));
    A = Wi.*P*P';
    B = Wi.*P;
    
    %matrix pT(x)
    x = points(i,1);
    y = points(i,2);
    Px = [1 x y x.^2 x.*y y.^2]';
    
    if det(A)>1e-30 
    %to avoid warning: matrix singular to machine precision  

        %shape function
        Nx(i,j) = Px'*(A\B);
  
        %with respect to x
        DWi1 = (ones(nP,1)*dw1(j));
        dAx1 = DWi1.*P*P';
        dAix1 = -A\dAx1/A;
        dBx1 = DWi1.*P;
        dPx1 = [0 1 0 2*x y 0]';
        dNx(i,j) = dPx1'*(A\B) + Px'*dAix1*B + Px'*(A\dBx1);
        
        %with respect to y
        DWi2 = (ones(nP,1)*dw2(j));
        dAx2 = DWi2.*P*P';
        dAix2 = -A\dAx2/A;
        dBx2 = DWi2.*P;
        dPx2 = [0 0 1 0 x 2*y]';
        dNy(i,j) = dPx2'*(A\B) + Px'*dAix2*B + Px'*(A\dBx2);
    
    else, det(A);
    end
end

end