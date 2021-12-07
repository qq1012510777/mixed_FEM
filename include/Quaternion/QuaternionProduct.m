function f = QuaternionProduct(Quaternion_t1, Quaternion_t2, Quaternion_t3, ...
        Quaternion_t4, ...
        Quaternion_f1, Quaternion_f2, Quaternion_f3, ...
        Quaternion_f4)

    t1 = GetVector(Quaternion_t1, Quaternion_t2, Quaternion_t3, ...
        Quaternion_t4);

    t2 = GetVector(Quaternion_f1, Quaternion_f2, Quaternion_f3, ...
        Quaternion_f4);

    scalar = Quaternion_t1 * Quaternion_f1 - dot(t1, t2);
    vector = Quaternion_t1 .* t2 + Quaternion_f1 .* t1 + cross(t1, t2);

    f = SetQuaternion(scalar, vector(1), vector(2), vector(3));

end
