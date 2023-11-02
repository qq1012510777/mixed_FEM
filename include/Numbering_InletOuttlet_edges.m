function f = Numbering_InletOuttlet_edges(dirichlet, nodes2edge)
    f = zeros(size(dirichlet, 1), 1);

    for i = 1:size(dirichlet, 1)
        node = dirichlet(i, [1, 2]);

        f(i) = nodes2edge(node(1), node(2));

    end

end
