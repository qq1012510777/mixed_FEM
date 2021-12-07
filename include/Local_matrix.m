function [alpha_ij, b_j] = Local_matrix(p1, p2, p3, global_vec, local_vec)
    xi = [0.0915762135,
        0.0915762135,
        0.8168475730,
        0.4459484909,
        0.4459484909,
        0.1081030182];
    eta = [0.8168475730,
        0.0915762135,
        0.0915762135,
        0.1081030182,
        0.4459484909,
        0.4459484909];
    w = [0.1099517437,
        0.1099517437,
        0.1099517437,
        0.2233815897,
        0.2233815897,
        0.2233815897];

    Area_ = Area_tri(p1, p2, p3);
    
    b_k = B_k(p1, p2, p3);

    phi = ["transpose((2^0.5) .* [x, y])";
        "transpose([-1 + x, y])";
        "transpose([x, y - 1])"];

    alpha_ij = zeros(3, 3);
    b_j = zeros(1, 3);

    for i = 1:6

        for j = 1:6
            xi_i = xi(i);
            eta_i = eta(i);

            xi_j = xi(j);
            eta_j = eta(j);

            w_i = w(i);
            w_j = w(j);

            %------------------alpha
            for ik = 1:3

                for jk = 1:3

                    edge_no_i = ik;
                    edge_no_j = jk;

                    phi_no_i = mod(ik + 1, 3) + 1;
                    phi_no_j = mod(jk + 1, 3) + 1;

                    x = xi_i;
                    y = eta_i;

                    eval(['phi_i = ', char(phi(phi_no_i)), ';']);

                    x = xi_j;
                    y = eta_j;

                    eval(['phi_j = ', char(phi(phi_no_j)), ';']);

                    sign_1 = 1;
                    sign_2 = 1;

                    global_vec_i = global_vec(edge_no_i, :);
                    local_vec_i = local_vec(edge_no_i, :);
                    tmp_1 = global_vec_i - local_vec_i;

                    if (norm(tmp_1) < 1e-4) % same direction
                        sign_1 = 1;
                    else
                        sign_1 = -1;
                    end

                    global_vec_j = global_vec(edge_no_j, :);
                    local_vec_j = local_vec(edge_no_j, :);
                    tmp_2 = global_vec_j - local_vec_j;

                    if (norm(tmp_2) < 1e-4)
                        sign_2 = -1;
                    else
                        sign_2 = 1;
                    end
                    
                    alpha_ij(ik, jk) = alpha_ij(ik, jk) + w_i * w_j * 1 / abs(det(b_k)) * dot(b_k * phi_i, b_k * phi_j);
                end

            end

            %-----------------------------
        end

    end

    %------------------b
    partial_phi = [2 * 2^0.5, 2, 2];

    for i = 1:3
        edge_no_i = i;

        phi_no_i = mod(ik + 1, 3) + 1;

        sign_1 = 1;
        global_vec_i = global_vec(edge_no_i, :);
        local_vec_i = local_vec(edge_no_i, :);
        tmp_1 = global_vec_i - local_vec_i;

        if (norm(tmp_1) < 1e-4) % same direction
            sign_1 = 1;
        else
            sign_1 = -1;
        end
        
        b_j(i) = (-1) / abs(det(b_k)) * sign_1 * partial_phi(phi_no_i) * Area_;
    end

end
