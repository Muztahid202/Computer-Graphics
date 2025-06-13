#include <iostream>
#include <fstream>
#include <vector>
#include <stack>
#include <cmath>
#include <iomanip>

using namespace std;

// Structure to represent a 3D point/vector
struct Point {
    double x, y, z, w;

    Point() : x(0), y(0), z(0), w(1) {}
    Point(double x, double y, double z) : x(x), y(y), z(z), w(1) {}
};

// Structure to represent a 4x4 matrix
struct Matrix {
    double m[4][4];

    Matrix() {
        // Initialize as identity matrix
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                m[i][j] = (i == j) ? 1.0 : 0.0;
            }
        }
    }

    // Matrix multiplication
    Matrix operator*(const Matrix& other) const {
        Matrix result;
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                result.m[i][j] = 0;
                for (int k = 0; k < 4; k++) {
                    result.m[i][j] += m[i][k] * other.m[k][j];
                }
            }
        }
        return result;
    }

    // Transform a point using this matrix
    Point transformPoint(const Point& p) const {
        Point result;
        result.x = m[0][0] * p.x + m[0][1] * p.y + m[0][2] * p.z + m[0][3] * p.w;
        result.y = m[1][0] * p.x + m[1][1] * p.y + m[1][2] * p.z + m[1][3] * p.w;
        result.z = m[2][0] * p.x + m[2][1] * p.y + m[2][2] * p.z + m[2][3] * p.w;
        result.w = m[3][0] * p.x + m[3][1] * p.y + m[3][2] * p.z + m[3][3] * p.w;

        // Normalize if w is not 1
        if (result.w != 1.0 && result.w != 0.0) {
            result.x /= result.w;
            result.y /= result.w;
            result.z /= result.w;
            result.w = 1.0;
        }

        return result;
    }
};

Matrix createTranslationMatrix(double tx, double ty, double tz) {
    Matrix mat;
    mat.m[0][3] = tx;
    mat.m[1][3] = ty;
    mat.m[2][3] = tz;
    return mat;
}

Matrix createScalingMatrix(double sx, double sy, double sz) {
    Matrix mat;
    mat.m[0][0] = sx;
    mat.m[1][1] = sy;
    mat.m[2][2] = sz;
    return mat;
}

Point normalizeVector(const Point& v) {
    double length = sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
    if (length == 0) return v;
    return Point(v.x / length, v.y / length, v.z / length);
}

double dotProduct(const Point& a, const Point& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

Point crossProduct(const Point& a, const Point& b) {
    return Point(
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x
    );
}

// Rodrigues' rotation formula
Point rotateVector(const Point& x, const Point& a, double angle) {
    angle = angle * M_PI / 180.0; // Convert to radians
    double cosTheta = cos(angle);
    double sinTheta = sin(angle);

    // cosθ * x
    Point term1 = Point(x.x * cosTheta, x.y * cosTheta, x.z * cosTheta);
    
    // (1 - cosθ)(a.x) * a
    double dot = dotProduct(a, x);
    Point term2 = Point(a.x * (1 - cosTheta) * dot, 
                        a.y * (1 - cosTheta) * dot, 
                        a.z * (1 - cosTheta) * dot);
    
    // sinθ * (a × x)
    Point cross = crossProduct(a, x);
    Point term3 = Point(cross.x * sinTheta, cross.y * sinTheta, cross.z * sinTheta);
    
    // R(x) = term1 + term2 + term3
    return Point(term1.x + term2.x + term3.x,
                term1.y + term2.y + term3.y,
                term1.z + term2.z + term3.z);
}

Matrix createRotationMatrix(double angle, double ax, double ay, double az) {
    Point a = normalizeVector(Point(ax, ay, az));
    
    // Apply Rodrigues' formula to basis vectors
    Point c1 = rotateVector(Point(1, 0, 0), a, angle);
    Point c2 = rotateVector(Point(0, 1, 0), a, angle);
    Point c3 = rotateVector(Point(0, 0, 1), a, angle);
    
    Matrix mat;
    mat.m[0][0] = c1.x; mat.m[0][1] = c2.x; mat.m[0][2] = c3.x;
    mat.m[1][0] = c1.y; mat.m[1][1] = c2.y; mat.m[1][2] = c3.y;
    mat.m[2][0] = c1.z; mat.m[2][1] = c2.z; mat.m[2][2] = c3.z;
    
    return mat;

    
}

Matrix createViewMatrix(double eyeX, double eyeY, double eyeZ,
                       double lookX, double lookY, double lookZ,
                       double upX, double upY, double upZ) {
    // Step 1: Calculate l, r, u vectors
    Point l(lookX - eyeX, lookY - eyeY, lookZ - eyeZ);
    l = normalizeVector(l);
    
    Point up(upX, upY, upZ);
    Point r = crossProduct(l, up);
    r = normalizeVector(r);
    
    Point u = crossProduct(r, l);
    
    // Step 2: Create translation matrix T
    Matrix T;
    T.m[0][3] = -eyeX;
    T.m[1][3] = -eyeY;
    T.m[2][3] = -eyeZ;
    
    // Step 3: Create rotation matrix R
    Matrix R;
    R.m[0][0] = r.x; R.m[0][1] = r.y; R.m[0][2] = r.z;
    R.m[1][0] = u.x; R.m[1][1] = u.y; R.m[1][2] = u.z;
    R.m[2][0] = -l.x; R.m[2][1] = -l.y; R.m[2][2] = -l.z;
    
    // Step 4: View matrix V = R * T
    return R * T;
}

Matrix createProjectionMatrix(double fovY, double aspectRatio, double near, double far) {
    Matrix P;
    
    // Compute fovX and necessary values
    double fovX = fovY * aspectRatio;
    double t = near * tan(fovY * M_PI / 360.0); // fovY/2 in radians
    double r = near * tan(fovX * M_PI / 360.0); // fovX/2 in radians
    
    // Set projection matrix values
    P.m[0][0] = near / r;
    P.m[1][1] = near / t;
    P.m[2][2] = -(far + near) / (far - near);
    P.m[2][3] = -(2.0 * far * near) / (far - near);
    P.m[3][2] = -1.0;
    P.m[3][3] = 0.0;
    
    return P;
}