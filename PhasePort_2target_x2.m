% stream plot of 2-target networks in paper, before + after affine map
close all; clear all;
traj_green = [6, 92, 74]/255; 
traj_linewidth = 2; 
ss_orange = [250, 140, 0]/255; 
ss_MarkerSize = 60; 

k = [1    1    1    1]; % rate constants

% phase portraits of the 2-targets network: 
density = 0.4;      % density displayed
hhstream = 0.01;    % stepsize in t
C = 1;              % scaling of dx/dt and dy/dt for resolution of plot
gridsize = [0,3];

tt = gridsize(1):hhstream:gridsize(2);
[XX,YY] = meshgrid(tt, tt); 

ddx = k(1) - k(2).*XX.^3 - k(3).*XX.^3.*YY.^2 + k(4).*YY.^2 ; 
ddy = k(1) + k(2).*XX.^3 - k(3).*XX.^3.*YY.^2 - k(4).*YY.^2 ;
u = C*ddx; v = C*ddy;

figure 
h = streamslice(XX,YY, u,v, density);
set(h,'Color', traj_green, 'LineWidth', traj_linewidth)
axis([gridsize, gridsize]);
set(gca,'XTick',[], 'YTick', [])
hold on;  

% steady state: 
plot(1,1, '.', 'MarkerSize', ss_MarkerSize, 'Color', ss_orange); hold off;



% phase portraits of the 2-targets network *AFTER* affine map: 
density = 0.4;      % density displayed
hhstream = 0.01;    % stepsize in t

affine_map = [4/3, 1/2; 1/3, 1];    % affine transformation
Ay1 = [0;0];                        % vertices after affine map
Ay2 = affine_map*[3;0]; 
Ay3 = affine_map*[3;2]; 
Ay4 = affine_map*[0;2]; 
Ay5 = affine_map*[1;1]; 
Ay6 = affine_map*[2;1]; 

tt = gridsize(1):hhstream:gridsize(2);
[XX,YY] = meshgrid(tt, tt); 

ddx = k(1)*(Ay5(1)-Ay1(1)) + k(2)*XX.^4.*YY*(Ay6(1)-Ay2(1)) + k(3)*XX.^5.*YY.^3*(Ay6(1)-Ay3(1)) + k(4)*XX.*YY.^2*(Ay5(1)-Ay4(1)); 
ddy = k(1)*(Ay5(2)-Ay1(2)) + k(2)*XX.^4.*YY*(Ay6(2)-Ay2(2)) + k(3)*XX.^5.*YY.^3*(Ay6(2)-Ay3(2)) + k(4)*XX.*YY.^2*(Ay5(2)-Ay4(2));
u = C*ddx; v = C*ddy;

figure
h = streamslice(XX,YY, u,v, density);
set(h,'Color', traj_green, 'LineWidth', traj_linewidth)
axis([gridsize, gridsize]);
set(gca,'XTick',[], 'YTick', [])
hold on;  

% steady state: 
plot(1,1, '.', 'MarkerSize', ss_MarkerSize, 'Color', ss_orange); hold off;

