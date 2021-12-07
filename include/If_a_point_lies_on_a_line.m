function f = If_a_point_lies_on_a_line(p1, p2, p)
    slope_ = (p1(2) - p2(2)) / (p1(1) - p2(1));

    b = p1(2) - slope_ * p1(1);

    y_prime = slope_ * p(1) + b;

    if (abs(y_prime - p(2)) < 1e-3)
        f = 1;
    else
        f = 0;
    end

end
