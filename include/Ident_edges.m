function [edges2nodes, tris2edges] = Ident_edges(t)
    [NUM_T, ~] = size(t);

    edges2nodes = [];
    tris2edges = [];

    for i = 1:NUM_T

        for j = 1:3
            node1 = t(i, j);
            node2 = t(i, mod(j, 3) + 1);

            smallnode = 0;
            largenode = 0;

            if (node1 < node2)
                smallnode = node1; largenode = node2;
            elseif (node1 > node2)
                smallnode = node2; largenode = node1;
            else
                disp('wrong edge');
                exit
            end

            edge_tmp = [smallnode, largenode];
            
            lia = is_member_row(edge_tmp, edges2nodes);         
           
            if (lia == 0)
                [ml, ~] = size(edges2nodes);
                edges2nodes(ml + 1, :) = edge_tmp;
                tris2edges(i, j) = ml + 1;
            else

                tris2edges(i, j) = lia;
            end

        end

    end

end
