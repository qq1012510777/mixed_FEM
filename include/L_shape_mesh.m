function [coordinate element dirichlet dirichlet_att Neumann] = L_shape_mesh(min_edge, max_edge)

    model = createpde;
    geometryFromEdges(model, @lshapeg);
    mesh = generateMesh(model, 'Hmax', max_edge, 'Hmin', min_edge, 'GeometricOrder', 'linear');

    coordinate = mesh.Nodes';
    element = mesh.Elements';

    d1 = 0;

    % dirchilet
    for i = 1:size(element, 1)

        for j = 1:3
            node1 = element(i, j);
            node2 = element(i, mod(j, 3) + 1);

            x1 = coordinate(node1, 1);
            x2 = coordinate(node2, 1);

            if (x1 == -1 && x2 == -1)
                d1 = d1 + 1;

                dirichlet(d1, :) = [node1, node2];

                dirichlet_att(d1, :) = "in";
            end

        end

    end

    for i = 1:size(element, 1)

        for j = 1:3
            node1 = element(i, j);
            node2 = element(i, mod(j, 3) + 1);

            x1 = coordinate(node1, 1);
            x2 = coordinate(node2, 1);

            if (x1 == 1 && x2 == 1)
                d1 = d1 + 1;

                dirichlet(d1, :) = [node1, node2];

                dirichlet_att(d1, :) = "out";
            end

        end

    end

    % Neumann
    d1 = 0;

    for i = 1:size(element, 1)

        for j = 1:3
            node1 = element(i, j);
            node2 = element(i, mod(j, 3) + 1);

            y1 = coordinate(node1, 2);
            y2 = coordinate(node2, 2);

            if (y1 == -1 && y2 == -1)
                d1 = d1 + 1;

                Neumann(d1, :) = [node1, node2];

            end

        end

    end

    for i = 1:size(element, 1)

        for j = 1:3
            node1 = element(i, j);
            node2 = element(i, mod(j, 3) + 1);

            y1 = coordinate(node1, 2);
            y2 = coordinate(node2, 2);

            x1 = coordinate(node1, 1);
            x2 = coordinate(node2, 1);

            if ((y1 == 1 && y2 == 1 && x1 >= 0 && x1 <= 1 && x2 >= 0 && x2 <= 1) || ...
                    (y1 == 0 && y2 == 0 && x1 >= -1 && x1 <= 0 && x2 >= -1 && x2 <= 0) || ...
                    (x1 == 0 && x2 == 0 && y1 >= 0 && y1 <= 1 && y2 >= 0 && y2 <= 1))
                d1 = d1 + 1;

                Neumann(d1, :) = [node1, node2];

            end

        end

    end

end
