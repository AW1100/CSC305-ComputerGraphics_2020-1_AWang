#include "lab.hpp"

int main()
{
    try
    {
        // clang-format off
        std::array<float, 216> vertices
        {
            // Vertices          Colours
            -1.0f,-1.0f,-1.0f,  1.0f, 0.0f, 0.0f,
            -1.0f,-1.0f, 1.0f,  1.0f, 0.0f, 0.0f,
            -1.0f, 1.0f, 1.0f,  1.0f, 0.0f, 0.0f,
            -1.0f,-1.0f,-1.0f,  1.0f, 0.0f, 0.0f,
            -1.0f, 1.0f, 1.0f,  1.0f, 0.0f, 0.0f,
            -1.0f, 1.0f,-1.0f,  1.0f, 0.0f, 0.0f,

            1.0f, 1.0f,-1.0f,   1.0f, 1.0f, 0.0f,
            -1.0f,-1.0f,-1.0f,  1.0f, 1.0f, 0.0f,
            -1.0f, 1.0f,-1.0f,  1.0f, 1.0f, 0.0f,
            1.0f, 1.0f,-1.0f,   1.0f, 1.0f, 0.0f,
            1.0f,-1.0f,-1.0f,   1.0f, 1.0f, 0.0f,
            -1.0f,-1.0f,-1.0f,  1.0f, 1.0f, 0.0f,

            1.0f, 1.0f, 1.0f,   0.0f, 0.0f, 1.0f,
            1.0f,-1.0f,-1.0f,   0.0f, 0.0f, 1.0f,
            1.0f, 1.0f,-1.0f,   0.0f, 0.0f, 1.0f,
            1.0f,-1.0f,-1.0f,   0.0f, 0.0f, 1.0f,
            1.0f, 1.0f, 1.0f,   0.0f, 0.0f, 1.0f,
            1.0f,-1.0f, 1.0f,   0.0f, 0.0f, 1.0f,

            -1.0f, 1.0f, 1.0f,  0.0f, 1.0f, 0.0f,
            -1.0f,-1.0f, 1.0f,  0.0f, 1.0f, 0.0f,
            1.0f,-1.0f, 1.0f,   0.0f, 1.0f, 0.0f,
            1.0f, 1.0f, 1.0f,   0.0f, 1.0f, 0.0f,
            -1.0f, 1.0f, 1.0f,  0.0f, 1.0f, 0.0f,
            1.0f,-1.0f, 1.0f,   0.0f, 1.0f, 0.0f,

            1.0f, 1.0f, 1.0f,   0.0f, 1.0f, 1.0f,
            1.0f, 1.0f,-1.0f,   0.0f, 1.0f, 1.0f,
            -1.0f, 1.0f,-1.0f,  0.0f, 1.0f, 1.0f,
            1.0f, 1.0f, 1.0f,   0.0f, 1.0f, 1.0f,
            -1.0f, 1.0f,-1.0f,  0.0f, 1.0f, 1.0f,
            -1.0f, 1.0f, 1.0f,  0.0f, 1.0f, 1.0f,

            1.0f,-1.0f, 1.0f,   1.0f, 1.0f, 1.0f,
            -1.0f,-1.0f,-1.0f,  1.0f, 1.0f, 1.0f,
            1.0f,-1.0f,-1.0f,   1.0f, 1.0f, 1.0f, 
            1.0f,-1.0f, 1.0f,   1.0f, 1.0f, 1.0f,
            -1.0f,-1.0f, 1.0f,  1.0f, 1.0f, 1.0f,
            -1.0f,-1.0f,-1.0f,  1.0f, 1.0f, 1.0f

        };
        // clang-format on

        Program prog{1280, 720, "CSC305 Lab 6"};
        Triangle cube{};

        cube.loadShaders();
        cube.loadDataToGPU(vertices);

        prog.run(cube);

        prog.freeGPUData();
        cube.freeGPUData();
    }
    catch (OpenGLError& err)
    {
        fmt::print("OpenGL Error:\n\t{}\n", err.what());
    }

    return 0;
}
