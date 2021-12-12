function f = Numbering_sep_edges(element, nodes2edge)

    f = zeros(3 * size(element, 1), 2);

    for i = 1:size(element, 1)

        k = 1;
        for j = [2, 3, 1]
            node1 = element(i, j);
            node2 = element(i, mod(j, 3) + 1);

            edgeNO = nodes2edge(node1, node2);

            f((i - 1) * 3 + k, 1) = edgeNO;
            f((i - 1) * 3 + k, 2) = i;
            k = k + 1;
        end

    end

end
