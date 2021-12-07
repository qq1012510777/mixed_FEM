function f = B_k(p1, p2, p3)
    x1 = p1(1);
    y1 = p1(2);

    x2 = p2(1);
    y2 = p2(2);

    x3 = p3(1);
    y3 = p3(2);

    f = [x2 - x1, x3 - x1;
        y2 - y1, y3 - y1];
end
