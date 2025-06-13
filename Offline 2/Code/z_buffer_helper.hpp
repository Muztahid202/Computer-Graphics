#include "bitmap_image.hpp"
#include "matrix_op.hpp"

class Triangle
{
public:
    Point points[3];
    int color[3]; // RGB color

    Triangle() {
        for (int i = 0; i < 3; i++) {
            points[i] = Point();
            color[i] = 0; // Default color black
        }
    }


};

class Cofactor
{
    public:
        double c03;
        double c13;
        double c23;
        double c33;
    
        Cofactor(double c03 = 0, double c13 = 0, double c23 = 0, double c33 = 0)
            : c03(c03), c13(c13), c23(c23), c33(c33) {}
};

struct IntersectionPoint
{
    double x1;
    double x2;
};

//storing the triangle objects from the stage3.txt file
vector<Triangle> readTrianglesFromFile(const string& filename)
{
    ifstream stage3File(filename);
    
    if (!stage3File.is_open())
    {
        cerr << "Error: Could not open stage3 file" << endl;
        exit(1);
    }
    

    vector<Triangle> triangles;
    Triangle triangle;
    while (stage3File >> triangle.points[0].x >> triangle.points[0].y >> triangle.points[0].z
                      >> triangle.points[1].x >> triangle.points[1].y >> triangle.points[1].z
                      >> triangle.points[2].x >> triangle.points[2].y >> triangle.points[2].z
                     )
    {
        // Assign random color to the triangle
        triangle.color[0] = rand() % 256; // Red
        triangle.color[1] = rand() % 256; // Green
        triangle.color[2] = rand() % 256; // Blue

        triangles.push_back(triangle);
    }

    stage3File.close();
    return triangles;
}

//read from the config file
void readConfig(const string& filename, int& screenWidth, int& screenHeight,
                double& leftLimitX, double& bottomLimitY, double& zFront, double& zRear)
{
    ifstream configFile(filename);
    if (!configFile.is_open())
    {
        cerr << "Error: Could not open config file" << endl;
        exit(1);
    }

    configFile >> screenWidth >> screenHeight;

    double temp;
    configFile >> temp;
    leftLimitX = temp;
    configFile >> temp;
    bottomLimitY = temp;

    configFile >> zFront >> zRear;

    configFile.close();
}

//creating the matrix form for creating the equation og the plane from the 3 points of the triangle
void createMatrixFromTriangle(const Triangle& triangle, Matrix& mat)
{
    // Create a matrix from the triangle points
    mat.m[0][0] = triangle.points[0].x;
    mat.m[1][0] = triangle.points[0].y;
    mat.m[2][0] = triangle.points[0].z;
    mat.m[0][1] = triangle.points[1].x;
    mat.m[1][1] = triangle.points[1].y;
    mat.m[2][1] = triangle.points[1].z;
    mat.m[0][2] = triangle.points[2].x;
    mat.m[1][2] = triangle.points[2].y;
    mat.m[2][2] = triangle.points[2].z;
    mat.m[0][3] = 1.0; 
    mat.m[1][3] = 1.0;
    mat.m[2][3] = 1.0; 
    mat.m[3][0] = 1.0; 
    mat.m[3][1] = 1.0;
    mat.m[3][2] = 1.0;
    mat.m[3][3] = 1.0;

}

double getCofactor(const Matrix& mat, int row, int col)
{
    double cofactor[3][3];
    int r = 0, c = 0;

    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            if (i != row && j != col)
            {
                cofactor[r][c++] = mat.m[i][j];
                if (c == 3)
                {
                    c = 0;
                    r++;
                }
            }
        }
    }

    // Calculate determinant of the 3x3 matrix
    double det = cofactor[0][0] * (cofactor[1][1] * cofactor[2][2] - cofactor[1][2] * cofactor[2][1]) -
                 cofactor[0][1] * (cofactor[1][0] * cofactor[2][2] - cofactor[1][2] * cofactor[2][0]) +
                 cofactor[0][2] * (cofactor[1][0] * cofactor[2][1] - cofactor[1][1] * cofactor[2][0]);

    // Apply the cofactor sign
    return ((row + col) % 2 == 0 ? 1 : -1) * det;
}

//calculate z for given x and y
double calculateZ(double c03, double c13, double c23, double c33, double x, double y)
{
    double z = (c33 + c03 * x + c13 * y) / (-c23);
    return z;
}

//find the intersection points of a line with the triangle given y
IntersectionPoint findIntersectionPoints(const Triangle& triangle, double y)
{
    IntersectionPoint intersection = {0, 0};

    const Point& p0 = triangle.points[0]; // lowest y
    const Point& p1 = triangle.points[1]; // middle y
    const Point& p2 = triangle.points[2]; // highest y

    // If y is outside the triangle's vertical bounds, return invalid result
    if (y < p0.y || y > p2.y)
        return intersection;

    double x1, x2;

    // Helper lambda to interpolate x from two points at given y
    auto interpolateX = [&](const Point& a, const Point& b) -> double {
        double t = (y - a.y) / (b.y - a.y);
        return a.x + t * (b.x - a.x);
    };

    if (y < p1.y) {
        // Lower half → intersect edges (p0, p1) and (p0, p2)
        x1 = interpolateX(p0, p1);
        x2 = interpolateX(p0, p2);
    } else if (y > p1.y) {
        // Upper half → intersect edges (p1, p2) and (p0, p2)
        x1 = interpolateX(p1, p2);
        x2 = interpolateX(p0, p2);
    } else {
        // y == p1.y → exactly horizontal at middle point
        x1 = p1.x;
        x2 = interpolateX(p0, p2);
    }

    // Ensure x1 is the smaller one (left to right)
    if (x1 > x2) swap(x1, x2);

    intersection.x1 = x1;
    intersection.x2 = x2;

    return intersection;
}

int mapToPixelRow(double y, double Top_Y, double dy)
{
    int row = static_cast<int>((Top_Y-y)/dy);
    return row;
}

int mapToPixelCol(double x, double Left_X, double dx)
{
    int col =  static_cast<int>((x - Left_X)/dx);
    return col;
}

double mapToY(int row_no, double Top_Y, double dy)
{
    double y = Top_Y - (row_no * dy);
    return y;
}

double mapToX(int col_no, double Left_X, double dx)
{
    double x = Left_X + (col_no * dx);
    return x;
}

void scanLineAlgorithm(const vector<Triangle>& triangles, double Top_Y, double Bottom_Y,
    double Left_X, double Right_X, double dy, double dx,const vector<Cofactor>& cofactors, double zFront, double zRear,
                      double** zBuffer, bitmap_image& image)
{
    int triangleCount = 0;

    for(const auto& triangle : triangles)
    {
        double yhigh = triangle.points[2].y;
        double ylow = triangle.points[0].y;

        int topScaleLine; //row no
        int bottomScaleLine;

        //clipping the portion over the top and bottom of the screen
        if(yhigh > Top_Y) {
            topScaleLine = mapToPixelRow(Top_Y, Top_Y, dy);
        }
        else {
            topScaleLine = mapToPixelRow(yhigh, Top_Y, dy);
        }
        if(ylow < Bottom_Y)
        {
            bottomScaleLine = mapToPixelRow(Bottom_Y, Top_Y, dy);
        }
        else{
            bottomScaleLine = mapToPixelRow(ylow, Top_Y, dy);
        }

        for(int row = topScaleLine; row <= bottomScaleLine; row++)
        {
            double y = mapToY(row, Top_Y, dy);
            IntersectionPoint intersection = findIntersectionPoints(triangle, y);

            //clip on x axis
            if(intersection.x1 < Left_X) {
                intersection.x1 = Left_X;
            }
            if(intersection.x2 > Right_X) {
                intersection.x2 = Right_X;
            }

            int leftCol = mapToPixelCol(intersection.x1, Left_X, dx);
            int rightCol = mapToPixelCol(intersection.x2, Left_X, dx);
            
            for (int col = leftCol; col <= rightCol; col++)
            {
                double x = mapToX(col, Left_X, dx);

                //calculate the z
                double z = calculateZ(cofactors[triangleCount].c03,
                                       cofactors[triangleCount].c13,
                                       cofactors[triangleCount].c23,
                                       cofactors[triangleCount].c33, x, y);
                if( z >= zFront && z <= zRear )
                {
                    //updating the zbuffer 
                      if (z < zBuffer[row][col]) {
                       // cout<<"entered zBuffer update"<<endl;
                        zBuffer[row][col] = z;
                        // Update pixel color
                        image.set_pixel(col, row, 
                                       triangle.color[0], 
                                       triangle.color[1], 
                                       triangle.color[2]);
                    }
                }

            }

        }

        triangleCount++;
    }

}


void saveZBuffer(const string& filename, double** zBuffer, 
                 int screenWidth, int screenHeight, double zMax) {
    ofstream zBufferFile(filename);
    if (!zBufferFile.is_open()) {
        cerr << "Error: Could not create Z-buffer file" << endl;
        return;
    }
    
    zBufferFile << fixed << setprecision(6);
    
    for (int i = 0; i < screenHeight; i++) {
        for (int j = 0; j < screenWidth; j++) {
            if (zBuffer[i][j] < zMax) {
                zBufferFile << zBuffer[i][j] << "\t";
            }
        }
        zBufferFile << endl;
    }
    
    zBufferFile.close();
}
