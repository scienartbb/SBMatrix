SBMatrix
=======

C++ Template based.
Dynamic Dimentional Matrix and Vector class

Mat22 multiply

    typedef SBMatrix<float, 2, 2> Mat22;
    Mat22 a = {1,2,0,1}, b={1,0,0,1};
    Mat22 c = a*b;
    printf(@"%s - %f %f %f %f", __FUNCTION__, c.m[0][0], c.m[0][1], c.m[1][0], c.m[1][1]);


Mat22, Vec2 multiply

    typedef SBMatrix<float, 2, 2> Mat22;
    typedef SBMatrix<float, 2, 1> Vec2;
    Mat22 a = {2,0,0,3};
    Vec2 b = {1,1};
    Vec2 c = a*b;
    printf(@"%s - %f %f ", __FUNCTION__, c.m[0][0], c.m[1][0]);



