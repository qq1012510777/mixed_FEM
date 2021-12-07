function f = is_member_row(A, B)
    [m, ~] = size(B);
    f = 0;

    for i = 1:m
       
        if (norm(A - B(i, :)) < 1e-4)
            f = i;
            break
        end

        if (norm(A - B(i, :)) > 1e-4 && i == m)
            f = 0;
            break
        end

    end

end
