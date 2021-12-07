function f = show_edges(edges2nodes, p)
    hold on
    [m, ~] = size(edges2nodes);

    for i = 1:m
        node1 = edges2nodes(i, 1);
        node2 = edges2nodes(i, 2);

        node = 0.5 .* ( p(node1, :) + p(node2, :));

        text(node(1), node(2), ['(', num2str(i), ')'], 'interpreter', 'latex', 'Color','blue');
        hold on
    end
    f = 0;
end