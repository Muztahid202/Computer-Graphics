// --- Includes ---
// Standard Headers
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

// OpenGL / GLUT Headers
#ifdef __APPLE__
#include <GLUT/glut.h> // Use GLUT framework on macOS
#else
#include <GL/glut.h> // Use standard GLUT location on Linux/Windows
#endif


// Camera position and orientation
GLfloat eyex = 8, eyey = 8, eyez = 8;          // Camera position coordinates
GLfloat centerx = 0, centery = 0, centerz = 0; // Look-at point coordinates
GLfloat upx = 0, upy = 1, upz = 0;             // Up vector coordinates

GLfloat MOVE_SPEED = 0.1;    // Speed of movement
GLfloat ROTATION_SPEED = 0.5; // in degrees

bool showVelocityArrow = true;
bool isSimulationRunning = true;
bool pressedReset = false; // Flag to check if reset key is pressed
float ballRadius = 0.1f;
float cubeSize = 4.0f; // Size of the cube (from -1 to 1 is 2 units)
float cubeHalfSize = cubeSize / 2.0f;

const float GRAVITY = -9.8f;            
const float RESTITUTION = 0.75f;        // Bounce coefficient (0.7-0.8)
const float VELOCITY_THRESHOLD = 0.01f; // Minimum velocity to consider ball at rest

// Ball properties
GLfloat ballPos[3] = {0.0f, -cubeHalfSize + ballRadius, 0.0f};
GLfloat ballVelocity[3] = {0.0f, 0.0f, 0.0f};
GLfloat ballRotation[3] = {0.0f, 0.0f, 0.0f}; // Rotation angles (degrees) for each axis
float initialSpeed = 2.0f;                    // Initial speed when reset

// --- Function Declarations ---
void initGL();
void display();
void reshapeListener(GLsizei width, GLsizei height);
void keyboardListener(unsigned char key, int x, int y);
void specialKeyListener(int key, int x, int y);
void drawCube();

/**
 * Initialize OpenGL settings
 * Sets up background color and enables depth testing
 */
void initGL()
{
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f); // Black background
    glEnable(GL_DEPTH_TEST);              // Enable depth testing for z-culling
}

void drawBall()
{
    int slices = 32; // Number of vertical slices 
    int stacks = 16; // Number of horizontal stacks 

    glPushMatrix();
    glTranslatef(ballPos[0], ballPos[1], ballPos[2]);

    // Apply rotation based on ball's rotation state
    glRotatef(ballRotation[0], 1.0f, 0.0f, 0.0f); // X-axis rotation
    glRotatef(ballRotation[1], 0.0f, 1.0f, 0.0f); // Y-axis rotation
    glRotatef(ballRotation[2], 0.0f, 0.0f, 1.0f); // Z-axis rotation

    for (int i = 0; i < slices; i++)
    {
        float theta1 = (float)i * 2.0f * M_PI / slices;
        float theta2 = (float)(i + 1) * 2.0f * M_PI / slices;

        // Alternate colors between slices
        if (i % 2 == 0)
            glColor3f(1.0f, 0.0f, 0.0f); // Red
        else
            glColor3f(0.0f, 0.0f, 1.0f); // Blue

        for (int j = 0; j < stacks; j++)
        {
            float phi1 = (float)j * M_PI / stacks;
            float phi2 = (float)(j + 1) * M_PI / stacks;

            float x1 = sinf(phi1) * cosf(theta1);
            float y1 = cosf(phi1);
            float z1 = sinf(phi1) * sinf(theta1);

            float x2 = sinf(phi2) * cosf(theta1);
            float y2 = cosf(phi2);
            float z2 = sinf(phi2) * sinf(theta1);

            float x3 = sinf(phi2) * cosf(theta2);
            float y3 = cosf(phi2);
            float z3 = sinf(phi2) * sinf(theta2);

            float x4 = sinf(phi1) * cosf(theta2);
            float y4 = cosf(phi1);
            float z4 = sinf(phi1) * sinf(theta2);

            glBegin(GL_QUADS);
            glVertex3f(x1 * ballRadius, y1 * ballRadius, z1 * ballRadius);
            glVertex3f(x2 * ballRadius, y2 * ballRadius, z2 * ballRadius);
            glVertex3f(x3 * ballRadius, y3 * ballRadius, z3 * ballRadius);
            glVertex3f(x4 * ballRadius, y4 * ballRadius, z4 * ballRadius);
            glEnd();
        }
    }

    glPopMatrix();
}

/**
 * Draw a checkered floor 
 */
void drawCheckeredFloor()
{
    int size = 20; // Number of tiles along one dimension
    float tileSize = cubeSize / size;
    bool colorToggle = true;

    glBegin(GL_QUADS);
    for (int x = -size / 2; x < size / 2; x++)
    {
        for (int z = -size / 2; z < size / 2; z++)
        {
            if (colorToggle)
            {
                glColor3f(1.0f, 1.0f, 1.0f); // White
            }
            else
            {
                glColor3f(0.0f, 0.0f, 0.0f); // Black
            }
            colorToggle = !colorToggle;

            float x1 = x * tileSize;
            float x2 = (x + 1) * tileSize;
            float z1 = z * tileSize;
            float z2 = (z + 1) * tileSize;

            glVertex3f(x1, -cubeHalfSize, z1);
            glVertex3f(x2, -cubeHalfSize, z1);
            glVertex3f(x2, -cubeHalfSize, z2);
            glVertex3f(x1, -cubeHalfSize, z2);
        }
        // Alternate starting color for each row
        if (size % 2 == 0)
            colorToggle = !colorToggle;
    }
    glEnd();
}

void drawVelocityArrow()
{
    if (!showVelocityArrow)
    {
        return;
    }

    glPushMatrix();
    glTranslatef(ballPos[0], ballPos[1], ballPos[2]);

    // Calculate direction (normalized) or use default if velocity is zero
    float dirX = 0.0f, dirY = 1.0f, dirZ = 0.0f; // Default upward direction
    float speed = sqrt(ballVelocity[0] * ballVelocity[0] +
                       ballVelocity[1] * ballVelocity[1] +
                       ballVelocity[2] * ballVelocity[2]);

    if (speed > 0.0001f) // Use actual direction if moving
    {
        dirX = ballVelocity[0] / speed;
        dirY = ballVelocity[1] / speed;
        dirZ = ballVelocity[2] / speed;
    }

    // Constant arrow dimensions
    const float arrowLength = 0.4f;                     // Total arrow length
    const float headLength = 0.15f;                     // Length of arrowhead
    const float shaftLength = arrowLength - headLength; // Length of shaft
    const float headRadius = 0.03f;                     // Radius of arrowhead base

    // Draw arrow shaft (extend slightly into the arrowhead for seamless connection)
    glLineWidth(1.8f);
    glColor3f(0.0f, 1.0f, 0.0f); // Green arrow
    glBegin(GL_LINES);
    glVertex3f(0, 0, 0);
    // Extend shaft by 10% into the head space for better visual connection
    glVertex3f(dirX * (shaftLength + headLength * 0.1),
               dirY * (shaftLength + headLength * 0.1),
               dirZ * (shaftLength + headLength * 0.1));
    glEnd();

    // Draw arrow head with perfect alignment
    glPushMatrix();
    // Position at the end of the shaft (minus the overlap we added)
    glTranslatef(dirX * shaftLength,
                 dirY * shaftLength,
                 dirZ * shaftLength);

    // Calculate proper rotation to align with velocity direction
    float pitch = -atan2(dirY, sqrt(dirX * dirX + dirZ * dirZ)) * 180.0f / M_PI;
    float yaw = atan2(dirX, dirZ) * 180.0f / M_PI;

    glRotatef(yaw, 0.0f, 1.0f, 0.0f);
    glRotatef(pitch, 1.0f, 0.0f, 0.0f);

    // Draw the cone pointing along positive Z-axis after rotations
    glutSolidCone(headRadius, headLength, 16, 2);
    glPopMatrix();

    glPopMatrix();
}


void drawCube()
{
    glBegin(GL_QUADS);

    // Top face (y = cubeHalfSize) - Green
    glColor3f(0.0f, 1.0f, 0.0f);
    glVertex3f(cubeHalfSize, cubeHalfSize, -cubeHalfSize);
    glVertex3f(-cubeHalfSize, cubeHalfSize, -cubeHalfSize);
    glVertex3f(-cubeHalfSize, cubeHalfSize, cubeHalfSize);
    glVertex3f(cubeHalfSize, cubeHalfSize, cubeHalfSize);

    // Front face (z = cubeHalfSize) - Red
    glColor3f(1.0f, 0.0f, 0.0f);
    glVertex3f(cubeHalfSize, cubeHalfSize, cubeHalfSize);
    glVertex3f(-cubeHalfSize, cubeHalfSize, cubeHalfSize);
    glVertex3f(-cubeHalfSize, -cubeHalfSize, cubeHalfSize);
    glVertex3f(cubeHalfSize, -cubeHalfSize, cubeHalfSize);

    // Back face (z = -cubeHalfSize) - Yellow
    glColor3f(1.0f, 1.0f, 0.0f);
    glVertex3f(cubeHalfSize, -cubeHalfSize, -cubeHalfSize);
    glVertex3f(-cubeHalfSize, -cubeHalfSize, -cubeHalfSize);
    glVertex3f(-cubeHalfSize, cubeHalfSize, -cubeHalfSize);
    glVertex3f(cubeHalfSize, cubeHalfSize, -cubeHalfSize);

    // Left face (x = -cubeHalfSize) - Blue
    glColor3f(0.0f, 0.0f, 1.0f);
    glVertex3f(-cubeHalfSize, cubeHalfSize, cubeHalfSize);
    glVertex3f(-cubeHalfSize, cubeHalfSize, -cubeHalfSize);
    glVertex3f(-cubeHalfSize, -cubeHalfSize, -cubeHalfSize);
    glVertex3f(-cubeHalfSize, -cubeHalfSize, cubeHalfSize);

    // Right face (x = cubeHalfSize) - Magenta
    glColor3f(1.0f, 0.0f, 1.0f);
    glVertex3f(cubeHalfSize, cubeHalfSize, -cubeHalfSize);
    glVertex3f(cubeHalfSize, cubeHalfSize, cubeHalfSize);
    glVertex3f(cubeHalfSize, -cubeHalfSize, cubeHalfSize);
    glVertex3f(cubeHalfSize, -cubeHalfSize, -cubeHalfSize);

    glEnd();

    // Draw the checkered floor separately
    drawCheckeredFloor();
}


void display()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    gluLookAt(eyex, eyey, eyez, centerx, centery, centerz, upx, upy, upz);

    drawCube(); 
    drawBall();
    drawVelocityArrow();

    glutSwapBuffers();
}

/**
 * Window reshape callback
 * Handles window resizing and maintains aspect ratio
 */
void reshapeListener(GLsizei width, GLsizei height)
{
    // Prevent division by zero
    if (height == 0)
        height = 1;

    // Calculate aspect ratio
    GLfloat aspect = (GLfloat)width / (GLfloat)height;

    // Set viewport to cover entire window
    glViewport(0, 0, width, height);

    // Set up perspective projection
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    // 45-degree field of view, aspect ratio, near and far clipping planes
    gluPerspective(45.0f, aspect, 0.1f, 100.0f);
}


// Rodrigues Formula for rotation(is applicable for rotating around any arbitrary axis)
void rotateVector(float &x, float &y, float &z, float ax, float ay, float az, float angle)
{
    // Rotate vector (x,y,z) about axis (ax,ay,az) by specified angle (degrees)
    float rad = angle * M_PI / 180.0f;
    float cosTheta = cos(rad);
    float sinTheta = sin(rad);

    // Normalize axis
    float len = sqrt(ax * ax + ay * ay + az * az);
    ax /= len;
    ay /= len;
    az /= len;

    // Rotation matrix components
    float m[3][3] = {
        {cosTheta + ax * ax * (1 - cosTheta), ax * ay * (1 - cosTheta) - az * sinTheta, ax * az * (1 - cosTheta) + ay * sinTheta},
        {ay * ax * (1 - cosTheta) + az * sinTheta, cosTheta + ay * ay * (1 - cosTheta), ay * az * (1 - cosTheta) - ax * sinTheta},
        {az * ax * (1 - cosTheta) - ay * sinTheta, az * ay * (1 - cosTheta) + ax * sinTheta, cosTheta + az * az * (1 - cosTheta)}};

    // Apply rotation
    float nx = m[0][0] * x + m[0][1] * y + m[0][2] * z;
    float ny = m[1][0] * x + m[1][1] * y + m[1][2] * z;
    float nz = m[2][0] * x + m[2][1] * y + m[2][2] * z;

    x = nx;
    y = ny;
    z = nz;
}

void updatePhysics(int value)
{
    if (isSimulationRunning)
    {
        float deltaTime = 0.016f; 

        // Apply gravity to vertical velocity
        ballVelocity[1] += GRAVITY * deltaTime;

        // Calculate displacement for this frame
        float dx = ballVelocity[0] * deltaTime;
        float dy = ballVelocity[1] * deltaTime;
        float dz = ballVelocity[2] * deltaTime;

        // Update position based on velocity
        ballPos[0] += dx;
        ballPos[1] += dy;
        ballPos[2] += dz;

        // Calculate rotation based on displacement (rolling effect)
        // Rotation angle in radians = distance traveled / radius
        float distance = sqrt(dx * dx + dz * dz);                // Horizontal distance only
        float angle = (distance / ballRadius) * (180.0f / M_PI); // Convert to degrees

        if (distance > 0.0001f)
        {
            // Calculate rotation axis (perpendicular to movement direction)
            float axisX = -dz / distance; // Rotate around Y axis for XZ movement
            float axisZ = dx / distance;

            // Update rotation (we'll store this as Euler angles for simplicity)
            ballRotation[0] += axisZ * angle; // Tilt based on movement direction
            ballRotation[2] += axisX * angle; // Roll based on movement direction

            // Keep angles within 0-360 range
            ballRotation[0] = fmod(ballRotation[0], 360.0f);
            ballRotation[1] = fmod(ballRotation[1], 360.0f);
            ballRotation[2] = fmod(ballRotation[2], 360.0f);
        }

        // Floor collision (y = -cubeHalfSize)
        if (ballPos[1] < (-cubeHalfSize + ballRadius))
        {
            ballPos[1] = -cubeHalfSize + ballRadius;          // Reposition to floor
            ballVelocity[1] = -ballVelocity[1] * RESTITUTION; // Bounce

            // Add some spin based on impact
            if (fabs(ballVelocity[0]) > 0.1f || fabs(ballVelocity[2]) > 0.1f)
            {
                ballRotation[0] += ballVelocity[2] * 10.0f; // Tilt based on Z velocity
                ballRotation[2] -= ballVelocity[0] * 10.0f; // Roll based on X velocity
            }

            // Stop the ball if velocity is very small
            if (fabs(ballVelocity[1]) < VELOCITY_THRESHOLD)
            {
                ballVelocity[1] = 0.0f;
            }
        }
        // Ceiling collision (y = cubeHalfSize)
        else if (ballPos[1] > (cubeHalfSize - ballRadius))
        {
            ballPos[1] = cubeHalfSize - ballRadius;           // Reposition to ceiling
            ballVelocity[1] = -ballVelocity[1] * RESTITUTION; // Bounce
        }

        // Left wall (x = -cubeHalfSize)
        if (ballPos[0] < (-cubeHalfSize + ballRadius))
        {
            ballPos[0] = -cubeHalfSize + ballRadius;
            ballVelocity[0] = -ballVelocity[0] * RESTITUTION;

            // Add spin on wall collision
            ballRotation[2] += ballVelocity[1] * 5.0f; // Spin around Z based on Y velocity
            ballRotation[1] += ballVelocity[2] * 5.0f; // Spin around Y based on Z velocity
        }
        // Right wall (x = cubeHalfSize)
        else if (ballPos[0] > (cubeHalfSize - ballRadius))
        {
            ballPos[0] = cubeHalfSize - ballRadius;
            ballVelocity[0] = -ballVelocity[0] * RESTITUTION;

            // Add spin on wall collision
            ballRotation[2] += ballVelocity[1] * 5.0f; // Spin around Z based on Y velocity
            ballRotation[1] += ballVelocity[2] * 5.0f; // Spin around Y based on Z velocity
        }

        // Front wall (z = cubeHalfSize)
        if (ballPos[2] > (cubeHalfSize - ballRadius))
        {
            ballPos[2] = cubeHalfSize - ballRadius;
            ballVelocity[2] = -ballVelocity[2] * RESTITUTION;

            // Add spin on wall collision
            ballRotation[0] += ballVelocity[1] * 5.0f; // Spin around X based on Y velocity
            ballRotation[1] += ballVelocity[0] * 5.0f; // Spin around Y based on X velocity
        }
        // Back wall (z = -cubeHalfSize)
        else if (ballPos[2] < (-cubeHalfSize + ballRadius))
        {
            ballPos[2] = -cubeHalfSize + ballRadius;
            ballVelocity[2] = -ballVelocity[2] * RESTITUTION;

            // Add spin on wall collision
            ballRotation[0] += ballVelocity[1] * 5.0f; // Spin around X based on Y velocity
            ballRotation[1] += ballVelocity[0] * 5.0f; // Spin around Y based on X velocity
        }

        glutPostRedisplay();
    }

    // Register the timer again (25ms ~ 40fps)
    glutTimerFunc(25, updatePhysics, 0);
}

void keyboardListener(unsigned char key, int x, int y)
{
    // Calculate and normalize view direction vector
    GLfloat fx = centerx - eyex;
    GLfloat fy = centery - eyey;
    GLfloat fz = centerz - eyez;
    float len = sqrt(fx * fx + fy * fy + fz * fz);
    fx /= len;
    fy /= len;
    fz /= len;

    // Calculate and normalize right vector
    GLfloat rx = fy * upz - fz * upy;
    GLfloat ry = fz * upx - fx * upz;
    GLfloat rz = fx * upy - fy * upx;
    len = sqrt(rx * rx + ry * ry + rz * rz);
    rx /= len;
    ry /= len;
    rz /= len;

    // Re-normalize up vector to ensure orthogonality
    upx = ry * fz - rz * fy;
    upy = rz * fx - rx * fz;
    upz = rx * fy - ry * fx;
    len = sqrt(upx * upx + upy * upy + upz * upz);
    upx /= len;
    upy /= len;
    upz /= len;

    switch (key)
    {
    // Yaw controls
    case '1': // Look left (Yaw)
        rotateVector(fx, fy, fz, upx, upy, upz, ROTATION_SPEED);
        centerx = eyex + fx;
        centery = eyey + fy;
        centerz = eyez + fz;
        break;

    case '2': // Look right (Yaw)
        rotateVector(fx, fy, fz, upx, upy, upz, -ROTATION_SPEED);
        centerx = eyex + fx;
        centery = eyey + fy;
        centerz = eyez + fz;
        break;

    // Pitch controls
    case '3': // Look up (Pitch)
        rotateVector(fx, fy, fz, rx, ry, rz, ROTATION_SPEED);
        rotateVector(upx, upy, upz, rx, ry, rz, ROTATION_SPEED);
        centerx = eyex + fx;
        centery = eyey + fy;
        centerz = eyez + fz;
        break;

    case '4': // Look down (Pitch)
        rotateVector(fx, fy, fz, rx, ry, rz, -ROTATION_SPEED);
        rotateVector(upx, upy, upz, rx, ry, rz, -ROTATION_SPEED);
        centerx = eyex + fx;
        centery = eyey + fy;
        centerz = eyez + fz;
        break;

    // Roll controls
    case '5': // Tilt clockwise (Roll)
        rotateVector(upx, upy, upz, fx, fy, fz, ROTATION_SPEED);
        // Recalculating right vector after roll
        rx = fy * upz - fz * upy;
        ry = fz * upx - fx * upz;
        rz = fx * upy - fy * upx;
        break;

    case '6': // Tilt counterclockwise (Roll)
        rotateVector(upx, upy, upz, fx, fy, fz, -ROTATION_SPEED);
        // Recalculating right vector after roll
        rx = fy * upz - fz * upy;
        ry = fz * upx - fx * upz;
        rz = fx * upy - fy * upx;
        break;

    case 'w': // Move upward without changing reference point
        printf("center point" " %f %f %f\n", centerx, centery, centerz);
        eyex += upx * MOVE_SPEED;
        eyey += upy * MOVE_SPEED;
        eyez += upz * MOVE_SPEED;
        break;

    case 's': // Move downward without changing reference point
        eyex -= upx * MOVE_SPEED;
        eyey -= upy * MOVE_SPEED;
        eyez -= upz * MOVE_SPEED;
        break;
    case ' ': // Toggle simulation
        isSimulationRunning = !isSimulationRunning;
        if (!isSimulationRunning)
        {
            pressedReset = false; // Reset the flag when simulation stops
        }
        break;

    case 'r': // Reset ball (only when simulation is off)
        if (!isSimulationRunning)
        {
            // Random position on the floor
            ballPos[0] = ((float)rand() / RAND_MAX - 0.5f) * (cubeHalfSize - ballRadius);
            ballPos[1] = -cubeHalfSize + ballRadius; // On the floor
            ballPos[2] = ((float)rand() / RAND_MAX - 0.5f) * (cubeHalfSize - ballRadius);

            // Random upward velocity with minimum upward component
            float angle = (float)rand() / RAND_MAX * M_PI_4;             // 0-45 degrees from vertical
            float horizontalAngle = (float)rand() / RAND_MAX * 2 * M_PI; // Random direction

            // Calculate velocity components
            ballVelocity[0] = sinf(angle) * cosf(horizontalAngle) * initialSpeed;
            ballVelocity[1] = cosf(angle) * initialSpeed; // Always positive (upward)
            ballVelocity[2] = sinf(angle) * sinf(horizontalAngle) * initialSpeed;

            // Ensure minimum upward velocity
            if (ballVelocity[1] < initialSpeed * 0.5f)
            {
                ballVelocity[1] = initialSpeed * 0.5f;
            }

            // Reset rotation
            ballRotation[0] = 0.0f;
            ballRotation[1] = 0.0f;
            ballRotation[2] = 0.0f;

            pressedReset = true; 
            printf("Ball reset. Current speed: %.2f\n", initialSpeed);
        }
        break;

    case '+': // Increase initial speed (only when simulation is off and after reset)
        if (!isSimulationRunning && pressedReset)
        {
            initialSpeed += 0.5f;
            if (initialSpeed > 20.0f)
                initialSpeed = 20.0f; // Upper limit
            printf("Speed increased to: %.2f\n", initialSpeed);

            // Update current velocity while keeping direction
            float speed = sqrt(ballVelocity[0] * ballVelocity[0] +
                               ballVelocity[1] * ballVelocity[1] +
                               ballVelocity[2] * ballVelocity[2]);
            if (speed > 0)
            {
                float scale = initialSpeed / speed;
                ballVelocity[0] *= scale;
                ballVelocity[1] *= scale;
                ballVelocity[2] *= scale;
            }
        }
        break;

    case '-': // Decrease initial speed (only when simulation is off and after reset)
        if (!isSimulationRunning && pressedReset)
        {
            initialSpeed -= 0.5f;
            if (initialSpeed < 0.5f)
                initialSpeed = 0.5f; // Lower limit
            printf("Speed decreased to: %.2f\n", initialSpeed);

            // Update current velocity while keeping direction
            float speed = sqrt(ballVelocity[0] * ballVelocity[0] +
                               ballVelocity[1] * ballVelocity[1] +
                               ballVelocity[2] * ballVelocity[2]);
            if (speed > 0)
            {
                float scale = initialSpeed / speed;
                ballVelocity[0] *= scale;
                ballVelocity[1] *= scale;
                ballVelocity[2] *= scale;
            }
        }
        break;

    case 'v': // Toggle velocity arrow
        showVelocityArrow = !showVelocityArrow;
        break;

    case 27:
        exit(0);
        break; // ESC key
    }

    glutPostRedisplay();
}

/**
 * Special key input handler (arrow keys, function keys)
 * Provides camera orbit functionality
 */
void specialKeyListener(int key, int x, int y)
{

    // Calculate forward vector
    GLfloat fx = centerx - eyex;
    GLfloat fy = centery - eyey;
    GLfloat fz = centerz - eyez;
    float len = sqrt(fx * fx + fy * fy + fz * fz);
    fx /= len;
    fy /= len;
    fz /= len; // Normalize

    // Calculate right vector (cross product of forward and up)
    GLfloat rx = fy * upz - fz * upy;
    GLfloat ry = fz * upx - fx * upz;
    GLfloat rz = fx * upy - fy * upx;
    len = sqrt(rx * rx + ry * ry + rz * rz);
    rx /= len;
    ry /= len;
    rz /= len; // Normalize

    switch (key)
    {
    // the eye and the center of the object both moves along the forward vector in the same direction by the same amount(zoom in)
    case GLUT_KEY_UP: // Move forward
        eyex += fx * MOVE_SPEED;
        eyey += fy * MOVE_SPEED;
        eyez += fz * MOVE_SPEED;

        centerx += fx * MOVE_SPEED;
        centery += fy * MOVE_SPEED;
        centerz += fz * MOVE_SPEED;
        break;
    // the eye and the center of the object both moves along the forward vector in the opposite direction by the same amount(zoom out)
    case GLUT_KEY_DOWN: // Move backward
        eyex -= fx * MOVE_SPEED;
        eyey -= fy * MOVE_SPEED;
        eyez -= fz * MOVE_SPEED;

        centerx -= fx * MOVE_SPEED;
        centery -= fy * MOVE_SPEED;
        centerz -= fz * MOVE_SPEED;
        break;
    // we need to move along the right vector. here eye and center both will be changed by the same amount. as we are going along the right vector in the opposite direction so we are negating the right vector
    case GLUT_KEY_LEFT: // Move left
        eyex -= rx * MOVE_SPEED;
        eyey -= ry * MOVE_SPEED;
        eyez -= rz * MOVE_SPEED;

        centerx -= rx * MOVE_SPEED;
        centery -= ry * MOVE_SPEED;
        centerz -= rz * MOVE_SPEED;
        break;
    // we need to move along the right vector. here eye and center both will be changed by the same amount. as we are going along the right vector in the same direction so we are not negating the right vector
    case GLUT_KEY_RIGHT: // Move right
        eyex += rx * MOVE_SPEED;
        eyey += ry * MOVE_SPEED;
        eyez += rz * MOVE_SPEED;

        centerx += rx * MOVE_SPEED;
        centery += ry * MOVE_SPEED;
        centerz += rz * MOVE_SPEED;
        break;
    // here we need to move along the up vector.
    case GLUT_KEY_PAGE_UP: // Move up
        eyex += upx * MOVE_SPEED;
        eyey += upy * MOVE_SPEED;
        eyez += upz * MOVE_SPEED;

        centerx += upx * MOVE_SPEED; // Move the center point too
        centery += upy * MOVE_SPEED;
        centerz += upz * MOVE_SPEED;
        break;

    case GLUT_KEY_PAGE_DOWN: // Move down
        eyex -= upx * MOVE_SPEED;
        eyey -= upy * MOVE_SPEED;
        eyez -= upz * MOVE_SPEED;

        centerx -= upx * MOVE_SPEED; // Move the center point too
        centery -= upy * MOVE_SPEED;
        centerz -= upz * MOVE_SPEED;
        break;
    }
    glutPostRedisplay();
}

/**
 * Draw coordinate axes
 * X axis: red, Y axis: green, Z axis: blue
 */
void drawAxes()
{
    glLineWidth(3); // Set line thickness

    glBegin(GL_LINES);

    // X axis (red)
    glColor3f(1, 0, 0);
    glVertex3f(0, 0, 0);
    glVertex3f(1, 0, 0);

    // Y axis (green)
    glColor3f(0, 1, 0);
    glVertex3f(0, 0, 0);
    glVertex3f(0, 1, 0);

    // Z axis (blue)
    glColor3f(0, 0, 1);
    glVertex3f(0, 0, 0);
    glVertex3f(0, 0, 1);

    glEnd();
}

/**
 * Main function: Program entry point
 */
int main(int argc, char **argv)
{
    // Initialize GLUT
    glutInit(&argc, argv);

    // Configure display mode and window
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);
    glutInitWindowSize(1024, 768);
    glutInitWindowPosition(50, 50);
    glutCreateWindow("OpenGL 3D Drawing");

    // Register callback functions
    glutDisplayFunc(display);
    glutReshapeFunc(reshapeListener);
    glutKeyboardFunc(keyboardListener);
    glutSpecialFunc(specialKeyListener);

    // Initialize OpenGL settings
    initGL();

    // Seed random number generator
    srand(time(NULL));

    // Initialize ball with random position on the floor
    ballPos[0] = ((float)rand() / RAND_MAX - 0.5f) * (cubeHalfSize - ballRadius);
    ballPos[1] = -cubeHalfSize + ballRadius; // On the floor
    ballPos[2] = ((float)rand() / RAND_MAX - 0.5f) * (cubeHalfSize - ballRadius);

    // Register the physics timer
    glutTimerFunc(25, updatePhysics, 0);

    // Enter the GLUT event loop
    glutMainLoop();

    return 0;
}