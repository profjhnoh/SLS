#ifndef PRECISION_H
#define PRECISION_H

#include <complex>

// 실수 정밀도 설정
// CMake 옵션으로 USE_FLOAT 또는 USE_DOUBLE 정의됨
// 기본값은 DOUBLE

#ifdef USE_FLOAT
    typedef float Real;
    typedef std::complex<Real> ComplexReal;
    #define REAL_FORMAT "%.6f"
    #define PRECISION_NAME "float"
    #define REAL(x) ((Real)(x))

    // Eigen 타입 suffix (float precision)
    #define MatrixXcReal MatrixXcf
    #define VectorXcReal VectorXcf
    #define RowVectorXcReal RowVectorXcf
    #define Matrix4cReal Matrix4cf
    #define Vector4cReal Vector4cf
    #define Matrix2cReal Matrix2cf
    #define Vector2cReal Vector2cf
    #define VectorXReal VectorXf
    #define MatrixXReal MatrixXf
    #define ArrayXReal ArrayXf
    #define RowVectorXReal RowVectorXf
#else
    typedef double Real;
    typedef std::complex<Real> ComplexReal;
    #define REAL_FORMAT "%.12f"
    #define PRECISION_NAME "double"
    #define REAL(x) ((Real)(x))

    // Eigen 타입 suffix (double precision)
    #define MatrixXcReal MatrixXcd
    #define VectorXcReal VectorXcd
    #define RowVectorXcReal RowVectorXcd
    #define Matrix4cReal Matrix4cd
    #define Vector4cReal Vector4cd
    #define Matrix2cReal Matrix2cd
    #define Vector2cReal Vector2cd
    #define VectorXReal VectorXd
    #define MatrixXReal MatrixXd
    #define ArrayXReal ArrayXd
    #define RowVectorXReal RowVectorXd
#endif

#endif // PRECISION_H
