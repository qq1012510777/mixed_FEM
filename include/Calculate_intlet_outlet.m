function [inlet outlet] = Calculate_intlet_outlet(coordinate, dirichlet, dirichlet_att, nodes2edge, x)
    inlet = 0;
    outlet = 0;

    Num_inlet_eles = 0; Num_outlet_eles = 0;

    for i = 1:size(dirichlet, 1)
        node1 = dirichlet(i, 1);
        node2 = dirichlet(i, 2);

        edgeNO = nodes2edge(node1, node2);

        if (edgeNO <= 0)
            error('error in Calculate_intlet_outlet.m');
        end

        len = norm(coordinate(node1, :) - coordinate(node2, :));

        if (dirichlet_att(i, 1) == "in")
            inlet = inlet + abs(x(edgeNO, 1)) * len;
            Num_inlet_eles = Num_inlet_eles + 1;
            %inlet
            %len

        elseif (dirichlet_att(i, 1) == "out")
            outlet = outlet + abs(x(edgeNO, 1)) * len;
            Num_outlet_eles = Num_outlet_eles + 1;
            %outlet
            %len

        end

    end

    Num_inlet_eles
    Num_outlet_eles
end
