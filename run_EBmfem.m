clc;
close all
clear all

currentPath = fileparts(mfilename('fullpath'));

addpath(genpath([currentPath, '/include']));
addpath(genpath([currentPath, '/EBmfem']));

% A. Program EBMFEM for 2D Raviart-Thomas mixed finite element method
%    based on the edge-oriented basis function
%
%    C.Bahriawati and C. Carstensen,08-08-03
%    File <EBmfem.m>
%
%    M-files you need to run
%       <stimaB.m>, <edge.m>, <f.m>, <u_D.m>, <g.m> (optional)
%
%    Data-files you need to run
%       <coordinate.dat>, <element.dat>,
%       <Dirichlet.dat>, <Neumann.dat> (optional)
%
%    This program and corresponding data-files use Example 9.1 in
%    "Three Matlab Implementations of the Lowest-Order Raviart-Thomas
%     MFEM with a Posteriori Error Control" by C.Bahriawati and C. Carstensen

% A.1. The main program
% load coordinate.dat;
% load element.dat;
% load dirichlet.dat;
% load Neumann.dat;
%
% dirichlet = [2 3; 5 6; 6 7;];
% dirichlet_att = ["out"; "in"; "in"];
%
% Neumann = [7 8; 8 1; 1 2; 3 4; 4 5];

%[coordinate element dirichlet dirichlet_att Neumann] = L_shape_mesh(0.5, 0.7);
%[coordinate element dirichlet dirichlet_att Neumann] = Rectangle_mesh(10, 10);
[coordinate element dirichlet dirichlet_att Neumann] = Trapezoidal_mesh(20, 20);

[nodes2element, nodes2edge, noedges, edge2element, interioredge] = edge(element, coordinate);

% A.2. EBmfem
%function u=EBmfem(element,coordinate,dirichlet,Neumann,nodes2element,...
%         nodes2edge,noedges,edge2element);

% Assemble matrices B and C
B = sparse(noedges, noedges);
C = sparse(noedges, size(element, 1));

figure(3)
Show_mesh(coordinate, element, 1)
Show_edge_NO(coordinate, nodes2edge)
hold on
title('node NO (red) and edge NO (black)')

for j = 1:size(element, 1)
    coord = coordinate(element(j, :), :)';
    I = diag(nodes2edge(element(j, [2 3 1]), element(j, [3 1 2])));
    signum = ones(1, 3);
    signum(find(j == edge2element(I, 4))) = -1;
    B(I, I) = B(I, I) + diag(signum) * stimaB(coord) * diag(signum);
    n = coord(:, [3, 1, 2]) - coord(:, [2, 3, 1]);
    C(I, j) = diag(signum) * [norm(n(:, 1)) norm(n(:, 2)) norm(n(:, 3))]';
end

% Global stiffness matrix A
A = sparse(noedges + size(element, 1), noedges + size(element, 1));
A = [B, C, ;
    C', sparse(size(C, 2), size(C, 2))];
% Volume force
b = sparse(noedges + size(element, 1), 1);

for j = 1:size(element, 1)
    %f(sum(coordinate(element(j, :), :)) / 3)

    b(noedges + j) = -det([1, 1, 1; coordinate(element(j, :), :)']) * ...
    f(sum(coordinate(element(j, :), :)) / 3) / 6;
end

% Dirichlet conditions
for k = 1:size(dirichlet, 1)

    if (dirichlet_att(k, 1) == "out")
        b(nodes2edge(dirichlet(k, 1), dirichlet(k, 2))) = norm(coordinate(dirichlet(k, 1), :) - ...
            coordinate(dirichlet(k, 2), :)) * (u_D(sum(coordinate(dirichlet(k, :), :)) / 2) + 0);
    else
        b(nodes2edge(dirichlet(k, 1), dirichlet(k, 2))) = norm(coordinate(dirichlet(k, 1), :) - ...
            coordinate(dirichlet(k, 2), :)) * (u_D(sum(coordinate(dirichlet(k, :), :)) / 2) + 100);
    end

end

% Neumann conditions
if ~isempty(Neumann)
    tmp = zeros(noedges + size(element, 1), 1);
    tmp(diag(nodes2edge(Neumann(:, 1), Neumann(:, 2)))) = ...
        ones(size(diag(nodes2edge(Neumann(:, 1), Neumann(:, 2))), 1), 1);
    FreeEdge = find(~tmp);
    x = zeros(noedges + size(element, 1), 1);
    CN = coordinate(Neumann(:, 2), :) - coordinate(Neumann(:, 1), :);

    for j = 1:size(Neumann, 1)
        x(nodes2edge(Neumann(j, 1), Neumann(j, 2))) = ...
            0; % g(sum(coordinate(Neumann(j, :), :)) / 2, CN(j, :) * [0, -1; 1, 0] / norm(CN(j, :)));
        % g(sum(coordinate(Neumann(j, :), :)) / 2, CN(j, :) * [0, -1; 1, 0] / norm(CN(j, :)))
    end

    b = b - A * x;
    x(FreeEdge) = A(FreeEdge, FreeEdge) \ b(FreeEdge);
else
    x = A \ b;
end

figure(1)
ShowDisplacement(element, coordinate, x);
colorbar
p = fluxEB(element, coordinate, x, noedges, nodes2edge, edge2element);
colorbar
close 1
figure(2)
ShowFlux(element, coordinate, p);
pEval = fluxEBEval(element, coordinate, x, nodes2edge, edge2element); hold on
figure(2)
subplot(2, 1, 1); view(2);
title('Velocity x')
colorbar
subplot(2, 1, 2); view(2);
title('Velocity y')
colorbar
eta_T = Aposteriori(element, coordinate, dirichlet, Neumann, x, pEval);

% ------------------------my function-----------------------------
% ------------------------my function-----------------------------
% ------------------------my function-----------------------------

Pressure = x(size(x, 1) - size(element, 1) + 1:end);

Center_ele = [];

for i = 1:size(element, 1)
    node1 = element(i, 1);
    node2 = element(i, 2);
    node3 = element(i, 3);

    Center_ele(i, :) = (coordinate(node1, :) + coordinate(node2, :) + ...
        coordinate(node3, :)) .* (1/3);
end

% figure(4)
% title('Pressure'); hold on
% Show_mesh(coordinate, element, 0)
% hold on
% [x1 y1] = meshgrid(-1:0.2:1, -1:0.2:1);
% %[x1 y1] = meshgrid(0:5:80, 0:5:50);
% z1 = griddata(Center_ele(:, 1), Center_ele(:, 2), Pressure, x1, y1, 'v4');
% C1 = contour(x1, y1, z1);
% colorbar

[inlet outlet] = Calculate_intlet_outlet(coordinate, dirichlet, dirichlet_att, nodes2edge, x);
inlet
outlet
