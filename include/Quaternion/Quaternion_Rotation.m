function f = Quaternion_Rotation(Angle_degree, axis_x, axis_y, axis_z, x, y, z)

    Quaternion_t = NormalizeRotation(axis_x, axis_y, axis_z, Angle_degree);
    f = Rotation(x, y, z, Quaternion_t(1), Quaternion_t(2), Quaternion_t(3), ...
        Quaternion_t(4));
end
