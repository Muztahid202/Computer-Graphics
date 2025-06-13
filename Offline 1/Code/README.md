# Computer Graphics OpenGL Setup

This repository contains code and resources for CSE410 Computer Graphics assignments using OpenGL.

## Prerequisites
- Ubuntu Linux
- g++ compiler
- OpenGL libraries: `libglu1-mesa-dev`, `freeglut3-dev`, `mesa-common-dev`

## Installation
Run the following commands to install OpenGL dependencies:

```sh
sudo apt-get update
sudo apt-get install libglu1-mesa-dev freeglut3-dev mesa-common-dev
```

## Building the Project
To compile one of the example programs (e.g., `ball_simulation.cpp` or `clock.cpp`), use:

```sh
g++ ball_simulation.cpp -o ball_simulation -lglut -lGLU -lGL
g++ clock.cpp -o clock -lglut -lGLU -lGL
```

## Running the Program
After compiling, run the executable:

```sh
./ball_simulation
# or
./clock
```

## Project Structure
- `Code/` — Source code files
- `Spec/` — Assignment specifications