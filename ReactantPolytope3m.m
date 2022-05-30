% Given mtx of source complexes, plot Reactant Polytope of CRN in 3D
%   via triangulation 
% Given mtx of product complexes, plot reaction vectors as well 
close all; clear all; 
src_complx_col = [0,  0,  255]/255; %[0.4310    0.8085    0.3465]; 
target_complx_col = [117, 115, 255]/255; 
rxn_col = [0,  0,  255]/255; %[0.9359    0.8986    0.1081]; 
Newton_polytope_col = [252, 128, 165]/255; 
complex_size = 100; 
rxn_vec_width = 2; 




% rxns: v1 <=> v2,  v3 <=> v4 <=> v5
v1 = [1;0;0]; 
v2 = [2;0;0];

v3 = [1;1;0];
v4 = [0;0;1];
v5 = [0;1;0];

Ysrc_mtx = [v1,v2,v3,v4,v5];

Rxn_src = [v1,v2,v3,v4,v4,v5];
Rxn_targ = [v2,v1,v4,v3,v3,v4];
Rxn_vect = Rxn_targ - Rxn_src; 


% plotting Newton/reactant polytope: 
DT = delaunayTriangulation(transpose(Ysrc_mtx));
[C,v] = convexHull(DT);
trisurf(C,DT.Points(:,1),DT.Points(:,2),DT.Points(:,3), ...
       'FaceColor', Newton_polytope_col, 'FaceAlpha', 0.15, ...
       'EdgeColor', "black", 'EdgeAlpha', 0.5); 
hold on; 
grid off; 
title('Newton polytope of CRN')
xlabel('x'); ylabel('y'); zlabel('z')

% plotting complexes
figure(1); 
scatter3(Ysrc_mtx(1,:),Ysrc_mtx(2,:),Ysrc_mtx(3,:), ...
    complex_size, src_complx_col,'filled'); 

for i=1:length(Rxn_src)
    quiver3(Rxn_src(1,i),Rxn_src(2,i),Rxn_src(3,i), ...
        Rxn_vect(1,i),Rxn_vect(2,i),Rxn_vect(3,i), ...
        'Color', rxn_col, 'LineWidth', rxn_vec_width)
end


