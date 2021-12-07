function f = Rotation(x, y, z, ...
        Quaternion_t1, Quaternion_t2, Quaternion_t3, Quaternion_t4)

    t1 = SetQuaternion(0.0, x, y, z);

    t2 = QuaternionProduct(Quaternion_t1, Quaternion_t2, Quaternion_t3, Quaternion_t4, ...
        t1(1), t1(2), t1(3), t1(4));

    t3 = Conjugate_Q(Quaternion_t1, Quaternion_t2, Quaternion_t3, Quaternion_t4);

    t1 = QuaternionProduct(t2(1), t2(2), t2(3), t2(4), ...
        t3(1), t3(2), t3(3), t3(4));

    f = GetVector(t1(1), t1(2), t1(3), t1(4));

end
