function f = Show_edge_NO(p, edges)
    hold on

    num = size(edges, 1);

    for i = 1:num - 1

        for j = i + 1:num
            edgeNO = edges(i, j);

            if (edgeNO ~= 0)
                center1 = 0.5 * (p(i, :) + p(j, :));

                text(center1(1, 1), center1(1, 2), ['\textcircled{', num2str(edgeNO), '}'], 'interpreter', 'latex'); hold on
            end

        end

    end

end
