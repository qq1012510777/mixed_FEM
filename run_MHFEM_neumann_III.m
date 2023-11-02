clc;
close all
clear all

conductivity = -5;
conductivity = conductivity^(-1);

currentPath = fileparts(mfilename('fullpath'));

addpath(genpath([currentPath, '/include']));
addpath(genpath([currentPath, '/EBmfem']));

% [coordinate element dirichlet dirichlet_att Neumann] = L_shape_mesh(0.1, 0.15);
% [coordinate element dirichlet dirichlet_att Neumann] = Rectangle_mesh(10, 20);
% [coordinate element dirichlet dirichlet_att Neumann] = Trapezoidal_mesh(10, 15);
[coordinate element dirichlet dirichlet_att Neumann] = Funnel_mesh(4, 7); 

figure(3)
Show_mesh(coordinate, element, 1)

[nodes2element, nodes2edge, noedges, edge2element, interioredge] = edge(element, coordinate);

noelements = size(element, 1);
nointerioredge = size(interioredge, 1);

B = sparse(noelements * 3, noelements * 3);
C = sparse(noelements * 3, noelements);
D = sparse(noelements * 3, noedges)';

Sep_edge_NO = Numbering_sep_edges(element, nodes2edge);
% Inter_edge_NO = Numbering_interior_edges(interioredge, edge2element);
b = sparse(noelements * 3 + noelements + noedges, 1);
figure(3)
hold on
Show_edge_NO(coordinate, nodes2edge)
hold on
title('node NO (red) and edge NO (black)')

% ------ calculate Inlet length and outlet length
L_in = 0; L_out = 0;
for i = 1:size(dirichlet_att, 1)
    edgeNO = nodes2edge(dirichlet(i, 1), dirichlet(i, 2));

    ele = edge2element(edgeNO, 3);

    node = edge2element(edgeNO, [1, 2]);
    len = norm(coordinate(node(1), :) - coordinate(node(2), :));

    if (dirichlet_att(i, 1) == "in")
        L_in = L_in + len;
    else
        L_out = L_out + len;
    end
end

for j = 1:noelements
    coord = coordinate(element(j, :), :)';

    %--------------------
    I_t = [(j - 1) * 3 + 1; (j - 1) * 3 + 2; (j - 1) * 3 + 3];
    B(I_t, I_t) = B(I_t, I_t) + conductivity * stimaB(coord); %diag(signum) * stimaB(coord) * diag(signum);
    n = coord(:, [2, 3, 1]) - coord(:, [3, 1, 2]);
    C(I_t, j) = C(I_t, j) + [norm(n(:, 1)) norm(n(:, 2)) norm(n(:, 3))]'; % diag(signum) * [norm(n(:, 1)) norm(n(:, 2)) norm(n(:, 3))]';

    glob_edgeNO = Sep_edge_NO(I_t, 1);
    local_edgeNO_interi = find(edge2element(glob_edgeNO, 4) ~= 0); % 1 2 3

    for k = 1:size(local_edgeNO_interi, 1)
        lo_in_ed = I_t(local_edgeNO_interi(k, 1), 1);

        % NO of interior edge
        global_interior_edgeNO = glob_edgeNO(local_edgeNO_interi(k, 1), 1);
        D(global_interior_edgeNO, lo_in_ed) = D(global_interior_edgeNO, lo_in_ed) - ...
            norm(n(:, local_edgeNO_interi(k, 1)));
    end
    
    % neumann and non-flux
    local_edgeNO_interi = find(edge2element(glob_edgeNO, 4) == 0); % 1 2 3
    for k = 1:size(local_edgeNO_interi, 1)
        Node_NO = edge2element(glob_edgeNO(local_edgeNO_interi(k, 1)), [1 2]);

        lo_in_ed = I_t(local_edgeNO_interi(k, 1), 1);
        % NO of inletoutlet edge
        global_inletoutlet_edgeNO = glob_edgeNO(local_edgeNO_interi(k, 1), 1);

        if (ismember(Node_NO, dirichlet, 'rows'))
            
            D(global_inletoutlet_edgeNO, lo_in_ed) = D(global_inletoutlet_edgeNO, lo_in_ed) - ...
                norm(n(:, local_edgeNO_interi(k, 1)));
            
            b_opop = 0;
            index_local = find(ismember(dirichlet, Node_NO, 'rows')==1);

            if (dirichlet_att(index_local) == "in")
                b_opop = -norm(n(:, local_edgeNO_interi(k, 1)))/L_in;
            else
                b_opop = norm(n(:, local_edgeNO_interi(k, 1)))/L_out;
            end
            b(noelements * 3 + noelements + global_inletoutlet_edgeNO, 1) = b(noelements * 3 + noelements + global_inletoutlet_edgeNO, 1) ...
                - b_opop * 1e-5;

        elseif (ismember(Node_NO, Neumann, 'rows'))
            D(global_inletoutlet_edgeNO, lo_in_ed) = D(global_inletoutlet_edgeNO, lo_in_ed) - ...
                norm(n(:, local_edgeNO_interi(k, 1)));
        else
            error('undefined boundary edge')
        end
    end
end

D = D';

K = [B, C, D;
    C', sparse(noelements, noelements), sparse(noelements, noedges);
    D', sparse(noedges, noelements), sparse(noedges, noedges)];

%-----------------------------------------------
% for i = 1:size(dirichlet_att, 1)
%     edgeNO = nodes2edge(dirichlet(i, 1), dirichlet(i, 2));
% 
%     ele = edge2element(edgeNO, 3);
% 
%     node = edge2element(edgeNO, [1, 2]);
%     len = norm(coordinate(node(1), :) - coordinate(node(2), :));
% 
%     I = diag(nodes2edge(element(ele, [2 3 1]), element(ele, [3 1 2])));
% 
%     edgeNO_k = find(I == edgeNO) + (ele - 1) * 3;
% 
%     if (dirichlet_att(i, 1) == "in")
%         b(edgeNO_k, 1) = 100 * len;
%     else
%         b(edgeNO_k, 1) = 0 * len;
%     end
% 
% end
% for i = 1:size(Neumann, 1)
%     edgeNO = nodes2edge(Neumann(i, 1), Neumann(i, 2));
%     ele = edge2element(edgeNO, 3);
%     I = diag(nodes2edge(element(ele, [2 3 1]), element(ele, [3 1 2])));
%     edgeNO_k = find(I == edgeNO) + (ele - 1) * 3;
% 
%     for j = 1:size(K, 1)
%         K(:, edgeNO_k) = 0;
%         K(edgeNO_k, :) = 0;
%     end
% 
%     b(edgeNO_k, 1) = 0;
%     K(edgeNO_k, edgeNO_k) = 1;
% end

A_ = K(1:noelements * 3, 1:noelements * 3);
B_ = K(1:noelements * 3, noelements * 3 + 1:noelements * 3 + noelements);
C_ = K(1:noelements * 3, noelements * 3 + noelements + 1:end);
g_ = b(1:noelements * 3, :);
f_ = b(noelements * 3 + 1:noelements * 3 + noelements, :);
s_ = b(noelements * 3 + noelements + 1:end, :);

B_ = B_';
C_ = C_';

D_ = C_ * inv(A_) * C_' - C_ * inv(A_) * B_' * inv((B_ * inv(A_) * B_')) * B_ * inv(A_) * C_';
r_ = C_ * inv(A_) * g_ - s_ + C_ * inv(A_) * B_' * inv((B_ * inv(A_) * B_')) * (f_ - B_ * inv(A_) * g_);

pressure_edge = (D_) \ r_;
% 
pressure_ele = inv((B_ * inv(A_) * B_')) * (B_ * inv(A_) * g_ - B_ * inv(A_) * C_' * pressure_edge - f_);
% 
velocity_sep_edge = inv(A_) * (g_ - B_' * pressure_ele - C_' * pressure_edge);
% 
% velocity_edge = zeros(noedges, 1);
% xx = K\b;
% velocity_sep_edge = xx(1:noelements * 3);
% pressure_ele = xx(noelements * 3 + 1: noelements * 3 + noelements);
% pressure_ele_edge = xx(noelements * 3 + noelements + 1: end);
velocity_edge = zeros(noedges, 1);

for i = 1:size(Sep_edge_NO, 1)
    edgeNO = Sep_edge_NO(i, 1);

    eleNO = Sep_edge_NO(i, 2);

    if (edge2element(edgeNO, 3) == eleNO)
        velocity_edge(edgeNO, 1) = velocity_sep_edge(i, 1);
    end

end

Center_ = [];

for i = 1:size(edge2element, 1)
    Center_(i, :) = 0.5 * (coordinate(edge2element(i, 1), :) + coordinate(edge2element(i, 2), :));
end

edge_normal = Global_normal(edge2element, coordinate);

figure(4); title('Global normal vectors', 'interpreter', 'latex'); hold on
P = patch('Vertices', coordinate, 'Faces', element, 'FaceVertexCData', zeros(size(coordinate, 1), 1), 'FaceColor', 'interp', 'EdgeAlpha', 0.9, 'facealpha', 0);
hold on
quiver(Center_(:, 1), Center_(:, 2), edge_normal(:, 1), edge_normal(:, 2));
hold on

flux_normal = [];

for i = 1:size(edge2element, 1)
    vec_ = edge_normal(i, :);
    if (velocity_edge(i, :) < 0); vec_ = -vec_; end;
    flux_normal(i, :) = vec_ .* abs(velocity_edge(i, :));
end

figure(5); title('Flux normal vectors', 'interpreter', 'latex'); hold on
P = patch('Vertices', coordinate, 'Faces', element, 'FaceVertexCData', zeros(size(coordinate, 1), 1), 'FaceColor', 'interp', 'EdgeAlpha', 0.9, 'facealpha', 0);
hold on
quiver(Center_(:, 1), Center_(:, 2), flux_normal(:, 1), flux_normal(:, 2), 0.8, 'r');
hold on

[inlet outlet] = Calculate_intlet_outlet(coordinate, dirichlet, dirichlet_att, nodes2edge, velocity_edge);
inlet
outlet
figure(6); title('Pressure field', 'interpreter', 'latex'); hold on
patch('Vertices', coordinate, 'Faces', element, 'FaceVertexCData', full(pressure_ele), 'FaceColor', 'flat', 'EdgeAlpha', 0.9, 'facealpha', 1);
colorbar