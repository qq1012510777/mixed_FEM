function [coordinate element dirichlet dirichlet_att Neumann] = Funnel_mesh(min_edge, max_edge)

    model = createpde;

    Trapezoidal_x = [0, 40, 80, 80, 40, 0];
    Trapezoidal_y = [-70, -10, -70, 120, 10, 120];

    R1 = [3, 6, Trapezoidal_x, Trapezoidal_y]';
    gm = [R1];
    sf = 'R1';
    ns = char('R1');
    ns = ns';
    g = decsg(gm, sf, ns);
    geometryFromEdges(model, g);

    mesh = generateMesh(model, 'Hmax', max_edge, 'Hmin', min_edge, 'GeometricOrder', 'linear');
    

    % pdemesh(mesh)

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

            if (x1 == Trapezoidal_x(1) && x2 == Trapezoidal_x(1))
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

            if (x1 == Trapezoidal_x(3) && x2 == Trapezoidal_x(3))
                d1 = d1 + 1;

                dirichlet(d1, :) = [node1, node2];

                dirichlet_att(d1, :) = "out";
            end

        end

    end

    % Neumann
    d1 = 0;

    p1 = [Trapezoidal_x(1), Trapezoidal_y(1)];
    p2 = [Trapezoidal_x(2), Trapezoidal_y(2)];

    for i = 1:size(element, 1)

        for j = 1:3
            node1 = element(i, j);
            node2 = element(i, mod(j, 3) + 1);

            p_u1 = coordinate(node1, :);
            p_u2 = coordinate(node2, :);

            f_1 = If_a_point_lies_on_a_line(p1, p2, p_u1);
            f_2 = If_a_point_lies_on_a_line(p1, p2, p_u2);

            if (f_1 == 1 && f_2 == 1)
                d1 = d1 + 1;

                Neumann(d1, :) = [node1, node2];

            end

        end

    end

    p1 = [Trapezoidal_x(2), Trapezoidal_y(2)];
    p2 = [Trapezoidal_x(3), Trapezoidal_y(3)];

    for i = 1:size(element, 1)

        for j = 1:3
            node1 = element(i, j);
            node2 = element(i, mod(j, 3) + 1);

            p_u1 = coordinate(node1, :);
            p_u2 = coordinate(node2, :);

            f_1 = If_a_point_lies_on_a_line(p1, p2, p_u1);
            f_2 = If_a_point_lies_on_a_line(p1, p2, p_u2);

            if (f_1 == 1 && f_2 == 1)
                d1 = d1 + 1;

                Neumann(d1, :) = [node1, node2];

            end

        end

    end

    p1 = [Trapezoidal_x(4), Trapezoidal_y(4)];
    p2 = [Trapezoidal_x(5), Trapezoidal_y(5)];

    for i = 1:size(element, 1)

        for j = 1:3
            node1 = element(i, j);
            node2 = element(i, mod(j, 3) + 1);

            p_u1 = coordinate(node1, :);
            p_u2 = coordinate(node2, :);

            f_1 = If_a_point_lies_on_a_line(p1, p2, p_u1);
            f_2 = If_a_point_lies_on_a_line(p1, p2, p_u2);

            if (f_1 == 1 && f_2 == 1)
                d1 = d1 + 1;

                Neumann(d1, :) = [node1, node2];

            end

        end

    end

    p1 = [Trapezoidal_x(5), Trapezoidal_y(5)];
    p2 = [Trapezoidal_x(6), Trapezoidal_y(6)];

    for i = 1:size(element, 1)

        for j = 1:3
            node1 = element(i, j);
            node2 = element(i, mod(j, 3) + 1);

            p_u1 = coordinate(node1, :);
            p_u2 = coordinate(node2, :);

            f_1 = If_a_point_lies_on_a_line(p1, p2, p_u1);
            f_2 = If_a_point_lies_on_a_line(p1, p2, p_u2);

            if (f_1 == 1 && f_2 == 1)
                d1 = d1 + 1;

                Neumann(d1, :) = [node1, node2];

            end

        end

    end

end
