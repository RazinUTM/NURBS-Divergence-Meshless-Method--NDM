% Title: Code for 1D & 2D Gauss Integration
% Created by: Malaysian Society for Numerical Methods (MSNM)
% Copyright © 2024 . All Rights Reserved


function [xg,wg,fg] = Gauss_Int(xi,eta,ng,ON)

% remove redundant points
xx = unique(xi);
yy = unique(eta);

% 1D Gauss points & weight in natural [-1 , 1]
% source: https://gubner.ece.wisc.edu/gaussquad.pdf (Example 15)
u = sqrt(1./(4-1./(1:ng-1).^2)); % upper diag.
[V,Lambda] = eig(diag(u,1)+diag(u,-1));
[xg,i] = sort(diag(Lambda));
Vtop = V(1,:);
Vtop = Vtop(i);
wg = 2*Vtop.^2;


if length(xx)==1 && length(yy)==1 %unable to generate

xg = [0 0];       %gauss points
wg = 0;           %weight
fg = 0;           %flag by knot spans

elseif length(xx)==1 || length(yy)==1

if length(xx)==1, xx = [xx xx]; end
if length(yy)==1, yy = [yy yy]; end

xG = []; wG = [];
for j = 1:length(yy)-1
    for i = 1:length(xx)-1
        x1 = xx(i);
        x2 = xx(i+1);
        y1 = yy(j);
        y2 = yy(j+1);
    
        for k = 1:ng %over x- or y- direction
            N1 = 1/2 - 1/2*xg(k);
            N2 = 1/2 + 1/2*xg(k);
            xk = N1*x1 + N2*x2;
                        
            N1 = 1/2 - 1/2*xg(k);
            N2 = 1/2 + 1/2*xg(k);
            yk = N1*y1 + N2*y2;

            %weight multiplied with detJ
            if (x2-x1)==0
                wk = wg(k)*(y2-y1)/2;
            elseif (y2-y1)==0
                wk = wg(k)*(x2-x1)/2;
            end
            
            xG = cat(1,xG,[xk yk wk i*j]);
        end
    end 
end
xg = xG(:,1:2);   %gauss points
wg = xG(:,3);     %weight
fg = xG(:,4);     %flag by knot spans


else

%2D gauss integration
xG = []; wG = [];
for j = 1:length(yy)-1
    for i = 1:length(xx)-1
        x1 = xx(i);
        x2 = xx(i+1);
        y1 = yy(j);
        y2 = yy(j+1);
    
        for l = 1:ng %over y-direction
            for k = 1:ng %over x-direction
                N1 = 1/2 - 1/2*xg(k);
                N2 = 1/2 + 1/2*xg(k);
                xk = N1*x1 + N2*x2;
                
                N1 = 1/2 - 1/2*xg(l);
                N2 = 1/2 + 1/2*xg(l);
                yk = N1*y1 + N2*y2;
                
                %weight multiplied with detJ
                wk = wg(k)*wg(l).*(x2-x1)*(y2-y1)/(2*2);
                
                xG = cat(1,xG,[xk yk wk j i]);
            end
        end
    end 
end
xg = xG(:,1:2);   %gauss points
wg = xG(:,3);     %weight
fg = xG(:,4:5);   %flag by knot spans

end
%--------------------------------------
% Plotting | 1=ON | 0=OFF |
%--------------------------------------
xx = unique(xi);
yy = unique(eta);
if length(xx)==1, xx = [-1 0 1]; end
if length(yy)==1, yy = [-1 0 1]; end

if exist('ON','var')
figure
hold on
axis equal
mesh(xx,yy,yy'*xx*0)
plot(xG(:,1),xG(:,2),'xk')
axis tight
hold off
end

end