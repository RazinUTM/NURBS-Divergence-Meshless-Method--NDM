% Title: 2D SF Code considering NURBS
% Created by: Malaysian Society for Numerical Methods (MSNM)
% Copyright © 2024 . All Rights Reserved

function [Ri,dRds,dRdn] = NURBS_SF(points,cpw,p,q,xi,eta)
%{
    NURBS_SF (NURBS Shape Function)
    - receive control points & weights (cpw) at physical space
      and NURBS parameters: (p),(q),(xi),(eta)
    - transform arbitary coordinates (points):s,n 
      from parametric/natural space to physical space
    - return shape function [Ri] and the derivatives [dRix],[dRiy]
      and the mapped coordinates [xyi] and the Jacobian matrix [J]
%}

% NURBS SF (Ri) & derivatives (dRds,dRdn) in parametric (a,b) [NxM]
[NiMj,dNiMj,NidMj] = Bspline(points,p,q,xi,eta);

w = cpw(:,3); % weights
Ri = zeros(size(NiMj));
dRds = zeros(size(NiMj));
dRdn = zeros(size(NiMj));
for i = 1:size(points,1)
    
    %Terms W, dWs &dWn
    W = (NiMj(i,:)*w);
    dWs = (dNiMj(i,:)*w);
    dWn = (NidMj(i,:)*w);
    
    Ri(i,:) = NiMj(i,:).*w' ./ (NiMj(i,:)*w); %shape function
    dRds(i,:) = w'.*(dNiMj(i,:).*W - NiMj(i,:).*dWs)...
                        ./(NiMj(i,:)*w)^2; %derivative over xi
    dRdn(i,:) = w'.*(NidMj(i,:).*W - NiMj(i,:).*dWn)...
                        ./(NiMj(i,:)*w)^2; %derivative over eta
end

end

function [Ns,dNs,dNn] = Bspline(points,p,q,xi,eta)

%B-spline surface 
%Shape function & 1st derivatives
Ns = []; dNs = []; dNn = [];
for i = 1:size(points,1)
    a = points(i,1);
    b = points(i,2);
    
    %get the basis function (BF)
    [Ni,dNi] = Basis(a,p,xi);
    [Mj,dMj] = Basis(b,q,eta);
    
    %select appropriate BF based on knot span xi
    Nia = Ni(:,sum(xi(1:end-p-1)<=a));
    dNia = dNi(:,sum(xi(1:end-p-1)<=a));
    
    %select appropriate BF based on knot span eta
    Mjb = Mj(:,sum(eta(1:end-q-1)<=b));
    dMjb = dMj(:,sum(eta(1:end-q-1)<=b));
    
    %formulate SF (Ni*Mj) & derivatives (dNi*Nj,Ni*dMj)
    Nsi = Nia(:,1)*Mjb(:,1)';
    dNsi = dNia(:,1)*Mjb(:,1)';
    dNni = Nia(:,1)*dMjb(:,1)';
    
    %put all of them in matrix for output
    Ns = cat(1,Ns, Nsi(:)' );
    dNs = cat(1,dNs, dNsi(:)' );
    dNn = cat(1,dNn, dNni(:)' );
    
end

end

function [N_all,dN_all] = Basis(x,P,knotv)
%to calculate the basis functions at node x

nk = length(knotv); %no.of knot vector
ncp = nk-P-1; %no.of control points

N = zeros(ncp+1,P+1,nk);
dN = zeros(ncp+1,P+1,nk);

%for P=0
N(P+1:ncp,1,P+1:ncp)=eye(ncp-P);

%for P=1,2,3...
for j = 2:P+1
    pp = j-1;
    i = 1:ncp;
    
    %numerators (in vector)
    num1 = x-knotv(i)';
    num2 = knotv(i+pp+1)'-x;
    
    %denominators (in vector)
    den1 = (knotv(i+pp)-knotv(i))';
    den2 = (knotv(i+pp+1)-knotv(i+1))';
    
    %terms 1 & 2 (in vector)
    term1 = num1./den1; term1(~(den1))=0;
    term2 = num2./den2; term2(~(den2))=0;
    
    %terms 1 & 2 (for derivative) (in vector)
    term1d = pp./den1; term1d(~(den1))=0;
    term2d = -pp./den2; term2d(~(den2))=0;
    
    %calculate according to the formula (in vector)
    N(i,j,:) = term1.*N(i,j-1,:) + term2.*N(i+1,j-1,:);
    dN(i,j,:) = term1d.*N(i,j-1,:) + term2d.*N(i+1,j-1,:);
    
end
N_all = squeeze(N(1:end-1,P+1,:));
dN_all = squeeze(dN(1:end-1,P+1,:));

end
