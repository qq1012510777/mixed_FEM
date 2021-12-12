function f = Numbering_interior_edges(interioredge, edge2element)
    f = zeros(size(interioredge, 1), 1);

    for i = 1:size(interioredge, 1)
        node = interioredge(i, [1, 2]);

        for j = 1:size(edge2element, 1)

            if (edge2element(j, 1) == node(1) && edge2element(j, 2) == node(2))
                f(i, 1) = j;
                break;
            end

        end

        if (f(i, 1) == 0)
            error('wrong');
        end

    end

end
