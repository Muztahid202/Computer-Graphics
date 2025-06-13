#include <GL/glut.h>
#include <ctime>
#include <cmath>
#include <chrono>

constexpr float PI = 3.14159265358979323846f;

class AnalogClock {
private:
    std::chrono::time_point<std::chrono::system_clock> startTime;
    time_t lastMinute = 0;
    bool firstMinute = true;

public:
    void init() {
        glClearColor(0.0f, 0.0f, 0.0f, 1.0f); // Black background
        gluOrtho2D(-250, 250, -250, 250); 
        
        // Get current system time
        time_t rawtime;
        time(&rawtime);
        struct tm *timeinfo = localtime(&rawtime);
        lastMinute = timeinfo->tm_min;
        
        // Initialize startTime to align with current seconds
        startTime = std::chrono::system_clock::now() - 
                   std::chrono::seconds(timeinfo->tm_sec) - 
                   std::chrono::milliseconds(timeinfo->tm_sec * 1000 % 1000);
    }

    void drawCircle(float radius) {
        glBegin(GL_LINE_LOOP);
        for (int i = 0; i < 360; ++i) {
            float angle = i * PI / 180;
            glVertex2f(radius * cos(angle), radius * sin(angle));
        }
        glEnd();
    }

    void drawHand(float length, float width, float angle, float r, float g, float b) {
        glLineWidth(width);
        glColor3f(r, g, b);
        glBegin(GL_LINES);
        glVertex2f(0.0f, 0.0f);
        glVertex2f(length * sin(angle), length * cos(angle));
        glEnd();
    }

    void drawClockFace() {
        // outer circle
        glColor3f(1.0f, 1.0f, 1.0f); // White color
        drawCircle(200.0f);

        // hour markers
        for (int i = 0; i < 12; ++i) {
            float angle = i * 30 * PI / 180;
            glBegin(GL_LINES);
            glVertex2f(180 * sin(angle), 180 * cos(angle));
            glVertex2f(190 * sin(angle), 190 * cos(angle));
            glEnd();
        }

        // minute markers
        glColor3f(0.8f, 0.8f, 0.8f); // Light gray
        for (int i = 0; i < 60; ++i) {
            if (i % 5 != 0) { // Skip positions where hour markers are
                float angle = i * 6 * PI / 180;
                glBegin(GL_LINES);
                glVertex2f(185 * sin(angle), 185 * cos(angle));
                glVertex2f(190 * sin(angle), 190 * cos(angle));
                glEnd();
            }
        }
    }

    void updateClock() {
        // current system time is being fetched here
        time_t rawtime;
        time(&rawtime);
        struct tm *timeinfo = localtime(&rawtime);

        // smoothly seconds are being calculated here using high-resolution clock
        auto now = std::chrono::system_clock::now();
        auto elapsed = now - startTime;
        auto total_millis = std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count();
        float seconds = fmod(total_millis / 1000.0f, 60.0f);

        // CWe are checking if the minute has changed
        if (timeinfo->tm_min != lastMinute) {
            lastMinute = timeinfo->tm_min;
            startTime = now - std::chrono::milliseconds((int)(seconds * 1000) % 1000);
            firstMinute = false;
        }

        // Calculating angles
        float secondAngle = seconds * 6.0f * PI / 180.0f;
        float minuteAngle = timeinfo->tm_min * 6.0f * PI / 180.0f;
        float hourAngle = (timeinfo->tm_hour % 12 + timeinfo->tm_min / 60.0f) * 30.0f * PI / 180.0f;

        //After updating

        // Draw hour hand (white)
        drawHand(80.0f, 6.0f, hourAngle, 1.0f, 1.0f, 1.0f);

        // Draw minute hand (cyan)
        drawHand(120.0f, 4.0f, minuteAngle, 0.0f, 1.0f, 1.0f);

        // Draw second hand (yellow)
        drawHand(150.0f, 2.0f, secondAngle, 1.0f, 1.0f, 0.0f);

        // Draw center pin (red)
        glPointSize(6.0f);
        glBegin(GL_POINTS);
        glColor3f(1.0f, 0.0f, 0.0f);
        glVertex2f(0.0f, 0.0f);
        glEnd();
    }

    void display() {
        glClear(GL_COLOR_BUFFER_BIT);
        drawClockFace();
        updateClock();
        glFlush();
        glutSwapBuffers();
    }

    static void timerCallback(int value) {
        glutPostRedisplay();
        glutTimerFunc(16, timerCallback, 0); 
    }
};

AnalogClock myClock;

void displayWrapper() {
    myClock.display();
}

int main(int argc, char** argv) {
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
    glutInitWindowSize(500, 500);
    glutInitWindowPosition(100, 100);
    glutCreateWindow("OpenGL Analog Clock");
    
    myClock.init();
    glutDisplayFunc(displayWrapper);
    glutTimerFunc(0, AnalogClock::timerCallback, 0);
    
    glutMainLoop();
    return 0;
}