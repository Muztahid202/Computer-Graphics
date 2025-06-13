#include "z_buffer_helper.hpp"
#include <filesystem>

int main()
{

    filesystem::create_directory("Output");

    double eyeX, eyeY, eyeZ;
    double lookX, lookY, lookZ;
    double upX, upY, upZ;
    double fovY, aspectRatio, near, far;

    const string sceneFilePath = "../Test Cases/3/scene.txt";
    const string configFilePath = "../Test Cases/3/config.txt";

    // Stage 1: Model Transformation
    {
        ifstream sceneFile(sceneFilePath);
        ofstream stage1File("Output/stage1.txt");

        // Set precision for output
        stage1File << fixed << setprecision(7);

        // Stack for transformation matrices
        stack<Matrix> matrixStack;
        matrixStack.push(Matrix()); // Push identity matrix initially

        string command;
        while (sceneFile >> command)
        {
            if (command == "triangle")
            {
                // Read three points
                Point points[3];
                for (int i = 0; i < 3; i++)
                {
                    sceneFile >> points[i].x >> points[i].y >> points[i].z;
                    points[i].w = 1;
                }

                // Transform and output the points
                for (int i = 0; i < 3; i++)
                {
                    Point transformed = matrixStack.top().transformPoint(points[i]);
                    stage1File << transformed.x << " " << transformed.y << " " << transformed.z << endl;
                }
                stage1File << endl;
            }
            else if (command == "translate")
            {
                double tx, ty, tz;
                sceneFile >> tx >> ty >> tz;
                Matrix translation = createTranslationMatrix(tx, ty, tz);
                Matrix newTop = matrixStack.top() * translation;
                matrixStack.pop();
                matrixStack.push(newTop);
            }
            else if (command == "scale")
            {
                double sx, sy, sz;
                sceneFile >> sx >> sy >> sz;
                Matrix scaling = createScalingMatrix(sx, sy, sz);
                Matrix newTop = matrixStack.top() * scaling;
                matrixStack.pop();
                matrixStack.push(newTop);
            }
            else if (command == "rotate")
            {
                double angle, ax, ay, az;
                sceneFile >> angle >> ax >> ay >> az;
                Matrix rotation = createRotationMatrix(angle, ax, ay, az);
                Matrix newTop = matrixStack.top() * rotation;
                matrixStack.pop();
                matrixStack.push(newTop);
            }
            else if (command == "push")
            {
                matrixStack.push(matrixStack.top()); // Duplicate the top matrix
            }
            else if (command == "pop")
            {
                if (matrixStack.size() > 1)
                {
                    matrixStack.pop();
                }
            }
            else if (command == "end")
            {
                break;
            }
        }

        sceneFile.close();
        stage1File.close();
    }

    // Stage 2: View Transformation
    {
        // Reopen scene.txt for reading view parameters
        ifstream sceneFile(sceneFilePath);
        ifstream stage1File("Output/stage1.txt");
        ofstream stage2File("Output/stage2.txt");

        // Read gluLookAt parameters from scene.txt
        sceneFile >> eyeX >> eyeY >> eyeZ;
        sceneFile >> lookX >> lookY >> lookZ;
        sceneFile >> upX >> upY >> upZ;
        sceneFile >> fovY >> aspectRatio >> near >> far;

        // Create view transformation matrix
        Matrix V = createViewMatrix(eyeX, eyeY, eyeZ,
                                    lookX, lookY, lookZ,
                                    upX, upY, upZ);

        // Set precision for output
        stage2File << fixed << setprecision(7);

        // Read points from stage1.txt and apply view transformation
        double x, y, z;
        int pointCount = 0;

        while (stage1File >> x >> y >> z)
        {
            Point p(x, y, z);
            Point transformed = V.transformPoint(p);
            stage2File << transformed.x << " " << transformed.y << " " << transformed.z << endl;

            pointCount++;
            if (pointCount % 3 == 0)
            {
                stage2File << endl; // Add empty line after each triangle
            }
        }

        sceneFile.close();
        stage1File.close();
        stage2File.close();
    }

    // Stage 3: Projection Transformation
    {
        // Reopen scene.txt for reading perspective parameters
        ifstream sceneFile(sceneFilePath);
        ifstream stage2File("Output/stage2.txt");
        ofstream stage3File("Output/stage3.txt");

        // Create projection matrix
        Matrix P = createProjectionMatrix(fovY, aspectRatio, near, far);

        // Set precision for output
        stage3File << fixed << setprecision(7);

        // Read points from stage2.txt and apply projection transformation
        double x, y, z;
        int pointCount = 0;

        while (stage2File >> x >> y >> z)
        {
            Point p(x, y, z);
            Point transformed = P.transformPoint(p);
            stage3File << transformed.x << " " << transformed.y << " " << transformed.z << endl;

            pointCount++;
            if (pointCount % 3 == 0)
            {
                stage3File << endl; // Add empty line after each triangle
            }
        }

        sceneFile.close();
        stage2File.close();
        stage3File.close();
    }

    // Stage 4:  Clipping & scan conversion using Z-buffer algorithm
    {
        int screenWidth, screenHeight;
        double leftLimitX, rightLimitX, bottomLimitY, topLimitY, zFront, zRear;
        readConfig(configFilePath, screenWidth, screenHeight, leftLimitX, bottomLimitY, zFront, zRear);

        rightLimitX = -leftLimitX;
        topLimitY = -bottomLimitY;

        // Initialize Z-buffer
        double **zBuffer = new double *[screenHeight];
        for (int i = 0; i < screenHeight; i++)
        {
            zBuffer[i] = new double[screenWidth];
            for (int j = 0; j < screenWidth; j++)
            {
                zBuffer[i][j] = zRear; // Initialize with maximum depth
            }
        }

        // Initialize image
        bitmap_image image(screenWidth, screenHeight);
        image.set_all_channels(0, 0, 0); // Black background

        // Read triangles from config file
        vector<Triangle> triangles = readTrianglesFromFile("Output/stage3.txt");

        // creating a vector that will store the Cofactor objects
        vector<Cofactor> cofactors;

        // for every triangle object in the vector we will do the following
        for (int i = 0; i < triangles.size(); i++)
        {
            Triangle &triangle = triangles[i];

            Matrix mat;

            createMatrixFromTriangle(triangle, mat);

            // caluclate the cofactor for (0,3), (1,3), (2,3), (3,3)
            double c03 = getCofactor(mat, 0, 3);
            double c13 = getCofactor(mat, 1, 3);
            double c23 = getCofactor(mat, 2, 3);
            double c33 = getCofactor(mat, 3, 3);

            // create a cofactor object and push it to the vector
            Cofactor cofactor(c03, c13, c23, c33);
            cofactors.push_back(cofactor);
        }

        for (Triangle &triangle : triangles)
        {
            // Sort the points of the triangle based on y-coordinate
            sort(triangle.points, triangle.points + 3, [](const Point &a, const Point &b)
                 { return a.y < b.y; });
        }

        double dx = (rightLimitX - leftLimitX) / screenWidth;
        double dy = (topLimitY - bottomLimitY) / screenHeight;
        double Top_Y = topLimitY - (dy / 2.0);
        double Bottom_Y = bottomLimitY + (dy / 2.0);
        double Left_X = leftLimitX + (dx / 2.0);
        double Right_X = rightLimitX - (dx / 2.0);

        // scanLineAlgorithm to fill the zBuffer and image
        scanLineAlgorithm(triangles, Top_Y, Bottom_Y, Left_X, Right_X, dy, dx, cofactors, zFront, zRear, zBuffer, image);

        // Save outputs
        image.save_image("Output/out.bmp");
        saveZBuffer("Output/z-buffer.txt", zBuffer, screenWidth, screenHeight, zRear);

        // Free memory
        for (int i = 0; i < screenHeight; i++)
        {
            delete[] zBuffer[i];
        }
        delete[] zBuffer;
    }

    return 0;
}
