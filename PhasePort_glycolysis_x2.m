% stream plot of Selkov model of glycolysis in paper, before + after affine map
close all; clear all;
traj_green = [6, 92, 74]/255; 
traj_linewidth = 2; 
ss_orange = [250, 140, 0]/255; 
ss_MarkerSize = 60; 
osc_linewidth = 7.5;


% phase portraits of the reversible Selkov model: 
density = 0.6;      % density displayed
hhstream = 0.01;    % stepsize in t
C = 1;              % scaling of dx/dt and dy/dt for resolution of plot
gridsize = [0,2.75];

tt = gridsize(1):hhstream:gridsize(2);
[XX,YY] = meshgrid(tt, tt); 

Source_Mtx = [0,0;
        0,1;
        0,1;
        1,0;
        1,0;
        0,0
        2,1;
        3,0];
Target_Cplx = [0,1;
        0,0;
        1,0;
        0,1;
        0,0;
        1,0
        3,0;
        2,1];
Rxn_Vect = Target_Cplx - Source_Mtx; 
k = [0.5; 0.1;0.01;0.01;1;0.1;1;0.1]; 

ddx = zeros(size(XX));
ddy = zeros(size(XX));
for i = 1:length(Rxn_Vect(:,1))
    ddx = ddx + k(i)*(XX.^(Source_Mtx(i,1))).*(YY.^(Source_Mtx(i,2))).*Rxn_Vect(i,1); 
    ddy = ddy + k(i)*(XX.^(Source_Mtx(i,1))).*(YY.^(Source_Mtx(i,2))).*Rxn_Vect(i,2); 
end
ddx = C*ddx; 
ddy = C*ddy; 
u = ddx; v = ddy;

figure
hold on;  plot(0.423,1.773, '.', 'MarkerSize', ss_MarkerSize, 'Color', ss_orange); hold off;  %s.s.
h_orbit = streamline(XX,YY, u,v, 1, 0.715);
set(h_orbit,'Color', ss_orange, 'LineWidth', osc_linewidth)

h = streamslice(XX,YY, u,v, density);
set(h,'Color', traj_green, 'LineWidth', traj_linewidth)
axis([gridsize, gridsize]);
set(gca,'XTick',[], 'YTick', [])



% phase portraits of the reversible Selkov model *AFTER* rotating by pi/3: 

angle_rotate = pi/3; 
A = [cos(angle_rotate), -sin(angle_rotate); 
        sin(angle_rotate), cos(angle_rotate)];
for i = 1:length(Rxn_Vect(:,1))
    Source_Mtx(i,:) = A*Source_Mtx(i,:)'; 
    Target_Cplx(i,:) = A*Target_Cplx(i,:)'; 
end
Rxn_Vect = Target_Cplx - Source_Mtx; 


ddx = zeros(size(XX));
ddy = zeros(size(XX));
for i = 1:length(Rxn_Vect(:,1))
    ddx = ddx + k(i)*(XX.^(Source_Mtx(i,1))).*(YY.^(Source_Mtx(i,2))).*Rxn_Vect(i,1); 
    ddy = ddy + k(i)*(XX.^(Source_Mtx(i,1))).*(YY.^(Source_Mtx(i,2))).*Rxn_Vect(i,2); 
end
ddx = C*ddx; 
ddy = C*ddy; 
u = ddx; v = ddy;

figure
h = streamslice(XX,YY, u,v, density);
set(h,'Color', traj_green, 'LineWidth', traj_linewidth)
axis([gridsize, gridsize]);
set(gca,'XTick',[], 'YTick', [])

% steady states: 
A_inv = inv(A);
    old_ss = [0.423,1.773]; 
    new_ss_x = old_ss(1)^(A_inv(1,1))*old_ss(2)^(A_inv(2,1)); 
    new_ss_y = old_ss(1)^(A_inv(1,2))*old_ss(2)^(A_inv(2,2)); 
hold on;  
plot(new_ss_x, new_ss_y, '.', 'MarkerSize', ss_MarkerSize, 'Color', ss_orange); hold off;  %s.s.
