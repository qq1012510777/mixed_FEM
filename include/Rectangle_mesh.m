function [coordinate element dirichlet dirichlet_att Neumann] = Rectangle_mesh(min_edge, max_edge)

    model = createpde;
    R1 = [3, 4, 0, 80, 80, 0, 0, 0, 50, 50]';
    gm = [R1];
    sf = 'R1';
    ns = char('R1');
    ns = ns';
    g = decsg(gm, sf, ns);
    geometryFromEdges(model, g);

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

            if (x1 == 0 && x2 == 0)
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

            if (x1 == 80 && x2 == 80)
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

            if (y1 == 0 && y2 == 0)
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

            if (y1 == 50 && y2 == 50)
                d1 = d1 + 1;

                Neumann(d1, :) = [node1, node2];

            end

        end

    end

end
