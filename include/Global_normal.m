function f = Global_normal(edge2element, coordinate)
    f = zeros(size(edge2element, 1), 2);

    for i = 1:size(edge2element, 1)
        node = edge2element(i, 1:2);
        ty = (coordinate(node(2), :) - coordinate(node(1), :));
        ty = ty ./ norm(ty);
        f(i, :) = [ty(2), -ty(1)];
    end

end
