function f = Show_mesh(p, t, if_show_node_tags)
    [m, ~] = size(p);
    patch('Vertices', p, 'Faces', t, 'FaceVertexCData', zeros(m, 1), 'FaceColor', 'interp', 'EdgeAlpha', 0.9, 'facealpha', 0);
    hold on

    if (if_show_node_tags == 1)
        ty = [1:1:m];
        ty = string(ty);
        text(p(:, 1)', p(:, 2)', ty', 'Color','red', 'interpreter', 'latex', 'fontsize', 15);
    end
end