R, (T, qw, qx, qy, z, r1, r2, r3) = polynomial_ring(QQ, ["T", "qw", "qx", "qy","z", "r1", "r2", "r3"])
sys=[T^2-3, 4*qw^4+12*qw^2*qx^2+9*qx^4-6*qx^2*qy^2+qy^4-4*qw*qy*z-8*qw^2-12*qx^2+4*
qy^2-r1^2+z^2+4, -6*T*qw^2*qx*qy-8*T*qx*qy^3+2*T*qw*qx*z+4*qw^4+3*qw^2*qx^2+9*
qw^2*qy^2+12*qx^2*qy^2+4*qy^4+8*T*qx*qy+2*qw*qy*z-8*qw^2-8*qy^2-r2^2+z^2+4, 6*T
*qw^2*qx*qy+8*T*qx*qy^3-2*T*qw*qx*z+4*qw^4+3*qw^2*qx^2+9*qw^2*qy^2+12*qx^2*qy^2
+4*qy^4-8*T*qx*qy+2*qw*qy*z-8*qw^2-8*qy^2-r3^2+z^2+4, qw^2+qx^2+qy^2-1]
