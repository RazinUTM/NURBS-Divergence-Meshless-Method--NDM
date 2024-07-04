% Title: NDM loop (control volume) numbering
% Created by: Malaysian Society for Numerical Methods (MSNM)
% Copyright © 2024 . All Rights Reserved

function [S,C,P,W] = NDM_num(a,b,m,n,ng,ON)
% UDF to calculates
% - side nodes numbering/connectivity [S]
% - side conditions whether internal/external [C]
% - coordinates of points for integration [P]
% - length/weights for numerical integration [W]
% from given size domain (a) x (b), the no.of loops
% and no.of Gauss point per each side (ng)

xx = linspace(0,1,m+1);
yy = linspace(0,1,n+1);

s0 = (2*m*ng + 2*n*ng); %all external sides

% Numbering for horizontal sides
H = zeros(n+1,m*ng);
for i = 2:n %internal sides
    s1 = (s0+1):(s0+m*ng);
    H(i,:) = s1 + (i-2)*m*ng;
end
H(1,:) = 1:m*ng; %bottom side
H(end,:) = (m*ng:-1:1)+(m+n)*ng; %top side

(H(end:-1:1,:));

s0 = max(max(H)); %continue numbering from H
% Numbering for vertical sides
V = zeros(n*ng,m+1);
for i = 2:m %internal sides
    s1 = (s0+1):(s0+n*ng);
    V(:,i) = s1 + (i-2)*n*ng;
end
V(:,end) = (1:n*ng)+m*ng;
V(:,1) = (n*ng:-1:1)+(2*m+n)*ng;

(V(end:-1:1,:));

% Combine all sides
H1 = H(1:n,:)';
S1 = reshape(H1(:),ng,[])';
V1 = V(:,2:m+1);
S2 = reshape(V1(:),ng,[])';
A = reshape(1:m*n,[],m)';
S2 = S2(A(:),:);
H2 = H(2:end,:)';
S3 = reshape(H2(:),ng,[])';
S3 = S3(:,end:-1:1);
V2 = V(:,1:m);
S4 = reshape(V2(:),ng,[])';
S4 = S4(A(:),end:-1:1);

S = [S1 S2 S3 S4];

% Location of each side (based on the 9 conditions)
% 0=boundary || 1=internal
C = S(:,1:ng:end) > (2*m+2*n)*ng;

% Integration points to evaluate
[xg,wg] = Gauss_Points_1D(ng,1);

%generate the location of points
hx1 = repmat(xx(1:end-1),[ng 1]);
hx2 = reshape(hx1+xg*diff(xx),1,[]);
hx = repmat(hx2,[n+1,1]);
vx = repmat(xx,[n*ng,1]);

vy1 = repmat(yy(1:end-1),[ng 1]);
vy2 = reshape(vy1+xg*diff(yy),1,[]);
vy = repmat(vy2',[1,m+1]);
hy = repmat(yy,[m*ng,1])';

% Integration weights (side-length) to evaluate
hx2 = reshape(wg*diff(xx),1,[]);
wx = repmat(hx2,[n+1,1]);

vy2 = reshape(wg*diff(yy),1,[]);
wy = repmat(vy2',[1,m+1]);

% sorting & finalize all values
xy = sortrows([ ...
    [H(:) hx(:) hy(:) wx(:) wx(:)*0];    %horizontal points & weights
    [V(:) vx(:) vy(:) wy(:)*0 wy(:)] ]); %vertical points & weights
P = xy(:,2:3);
P = [P(:,1)*a , P(:,2)*b];
W = xy(:,4:5);
W = [W(:,1)*a , W(:,2)*b];
W = sum(W,2);


%--------------------------------------
% Plotting
%--------------------------------------
if exist('ON','var')
    h = figure(ON);
    set(h, 'WindowStyle', 'Docked');
    axis equal tight
    axis([0,a,0,b])
    hold on
    for i = 1:size(S,1) %elements (sub-area)
        x = P(S(i,end),1)+a/(2*m);
        y = P(S(i,1),2)+b/(2*n);
        x1 = min(P(S(i,:),1));
        x2 = max(P(S(i,:),1));
        y1 = min(P(S(i,:),2));
        y2 = max(P(S(i,:),2));
        X = [x1 x2 x2 x1 x1];
        Y = [y1 y1 y2 y2 y1];
        plot(X,Y,'c:')
        plot(P(S(i,:),1),P(S(i,:),2),'c:','LineWidth',1)
        text(x,y,num2str(i),'Color',[0, 0 ,1], ...
            "horizontalalignment", "center")
        pause(0.2)
    end
    for i = 1:size(P,1) %all integration points
        x = P(i,1);
        y = P(i,2);
        plot(x,y,'bx','markerSize',3,'LineWidth',1)
        text(x,y,num2str(i),'Color',[1, 0 ,0.5], ...
            "horizontalalignment", "center")
        pause(0.1)
    end
    hold off
end

end

function [xg,wg] = Gauss_Points_1D(ng,L)
% source: https://gubner.ece.wisc.edu/gaussquad.pdf (Example 15)
% This function
% calculates Gpoints [xg] & weight [wg] for 1D
% for a given no. of Gpoints ng x ng (ng)
% and the line length (L)

c = (1:ng-1)./sqrt(4*(1:ng-1).^2-1);
J = diag(c,-1) + diag(c,1); % diagonal values is 0
[V,D] = eig(J);
[X,ix] = sort(diag(D)); % sort diagonal eigenvalues

V(:,1:ng) = V(:,ix);
w = 2*V(1,:)'.^2;

xg = X*(L/2) + (L/2); % change of interval
wg = w*L/2;

end

