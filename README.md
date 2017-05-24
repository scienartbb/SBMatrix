SBMatrix
=======

C++ Template based.

Platform Independent. (iOS/Android/Mac/Win32/Linux)

N-Dim Matrix and Vector class

Simple, only 1 Header File.






Mat22 multiply

    typedef SBMatrix<float, 2, 2> Mat22;
    Mat22 a = {1,2,0,1}, b={1,0,0,1};
    Mat22 c = a*b;
    printf("%s - %f %f %f %f", __FUNCTION__, c.m[0][0], c.m[0][1], c.m[1][0], c.m[1][1]);


Mat22, Vec2 multiply

    typedef SBMatrix<float, 2, 2> Mat22;
    typedef SBMatrix<float, 2, 1> Vec2;
    Mat22 a = {2,0,0,3};
    Vec2 b = {1,1};
    Vec2 c = a*b;
    printf("%s - %f %f ", __FUNCTION__, c.m[0][0], c.m[1][0]);

You can use it Dynamic Dimentional.

    typedef SBMatrix<float, 3, 3> Mat3; // for 3x3
    typedef SBMatrix<float, 4, 4> Mat4; // for 4x4
    typedef SBMatrix<float, 5, 5> Mat5; // for 5x5

also non square!

    typedef SBMatrix<float, 3, 2> Mat32; // for 3x3
    typedef SBMatrix<float, 4, 3> Mat43; // for 4x4
    typedef SBMatrix<float, 5, 2> Mat52; // for 5x5

Matrix Inverse (Square Matrix)

    typedef SBMatrix<float, 2, 2> Mat22;
    Mat22 a = {2,0,0,3};
    Mat22 b = a.inv();

Others
- Determinent (det())
- Gauss Elimination (gausselim())
- Gauss Jordan Elimination (gjordelim())
- LU factorization (lu())
- QR factorization (qr())

Global Routines
- rotation
- translation
- scaling
- skewing
- lookAt
- projection (Orthogonal, Perspective)
- euler
- affine transform 2D / 3D
- lerp
- de casteljau
- bezier quad / cubic
- arc
- circle
- hermite
- perpendicular CW / CCW
- find unit circle
- find ellipse

