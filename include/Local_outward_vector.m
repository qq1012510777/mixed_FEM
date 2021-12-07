function [x y fx fy] = Local_outward_vector(p, t)
    [m, ~] = size(t);
    x = zeros(m, 3);
    y = zeros(m, 3);

    fx = zeros(m, 3);
    fy = zeros(m, 3);

    for i = 1:m

        for j = 1:3
            node1 = t(i, j);
            node2 = t(i, mod(j, 3) + 1);

            center = (p(node1, :) + p(node2, :)) * 0.5;
            
            x(i, j) = center(1);
            y(i, j) = center(2);
            
            direc_vec = p(node2, :) - p(node1, :);
            inward_vec = Quaternion_Rotation(90, 0, 0, 1, direc_vec(1), direc_vec(2), 0);
            outward_vec = -inward_vec;
            outward_vec = outward_vec ./ norm(outward_vec);
            fx(i, j) = outward_vec(1);
            fy(i, j) = outward_vec(2);
        end

    end

end
