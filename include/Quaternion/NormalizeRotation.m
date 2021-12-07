function f = NormalizeRotation(axis_x, axis_y, axis_z, angle_degree)
    axis_ = [axis_x, axis_y, axis_z];

    if (norm(axis_) < 1e-12)
        throw('Quaternion: the norm of the axis is too small, please chose a different one\n');
    end

    A = axis_ ./ norm(axis_);

    angle_radian = angle_degree * pi / 180.0;

    C = [cos(angle_radian / 2.0), ...
            A(1) * sin(angle_radian / 2.0), ...
            A(2) * sin(angle_radian / 2.0), ...
            A(3) * sin(angle_radian / 2.0)];

    f = C;
end
