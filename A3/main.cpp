#include "assignment.hpp"

static std::vector<float> generateSphereVertices()
{
    std::vector<float> vertices{};

    float radius = 0.50f;
    int sectorCount = 36;
    int stackCount = 18;
    float PI = 3.14159265f;
    glm::vec3 colour{ 1.0f,0.0f,0.0f };

    float x, y, z, xy;
    float nx, ny, nz, lengthInv = 1.0f / radius;

    float sectorStep = 2 * PI / sectorCount;
    float stackStep = PI / stackCount;
    float sectorAngle, stackAngle;

    for (int i = 0; i <= stackCount; ++i)
    {
        stackAngle = PI / 2 - i * stackStep;        // starting from pi/2 to -pi/2
        xy = radius * cosf(stackAngle);             // r * cos(u)
        z = radius * sinf(stackAngle);              // r * sin(u)

        // add (sectorCount+1) vertices per stack
        // the first and last vertices have same position and normal, but different tex coords
        for (int j = 0; j <= sectorCount; ++j)
        {
            sectorAngle = j * sectorStep;           // starting from 0 to 2pi

            // vertex position (x, y, z)
            x = xy * cosf(sectorAngle);             // r * cos(u) * cos(v)
            y = xy * sinf(sectorAngle);             // r * cos(u) * sin(v)
            vertices.push_back(x);
            vertices.push_back(y);
            vertices.push_back(z);

            vertices.push_back(colour.x);
            vertices.push_back(colour.y);
            vertices.push_back(colour.z);

            // normalized vertex normal (nx, ny, nz)
            nx = x * lengthInv;
            ny = y * lengthInv;
            nz = z * lengthInv;
            vertices.push_back(nx);
            vertices.push_back(ny);
            vertices.push_back(nz);

        }
    }

    return vertices;
}


int main()
{
    try
    {
        // clang-format off
        std::vector<float> vertices
        {
            // Vertices          Colours                normal
            -0.5f, -0.5f, -0.5f,  1.0f, 0.0f, 0.0f,   0.0f,  0.0f, -1.0f,
             0.5f, -0.5f, -0.5f,  1.0f, 0.0f, 0.0f,   0.0f,  0.0f, -1.0f,
             0.5f,  0.5f, -0.5f,  1.0f, 0.0f, 0.0f,   0.0f,  0.0f, -1.0f,
             0.5f,  0.5f, -0.5f,  1.0f, 0.0f, 0.0f,   0.0f,  0.0f, -1.0f,
            -0.5f,  0.5f, -0.5f,  1.0f, 0.0f, 0.0f,   0.0f,  0.0f, -1.0f,
            -0.5f, -0.5f, -0.5f,  1.0f, 0.0f, 0.0f,   0.0f,  0.0f, -1.0f,
                                                     
            -0.5f, -0.5f,  0.5f,  1.0f, 1.0f, 0.0f,   0.0f,  0.0f,  1.0f,
             0.5f, -0.5f,  0.5f,  1.0f, 1.0f, 0.0f,   0.0f,  0.0f,  1.0f,
             0.5f,  0.5f,  0.5f,  1.0f, 1.0f, 0.0f,   0.0f,  0.0f,  1.0f,
             0.5f,  0.5f,  0.5f,  1.0f, 1.0f, 0.0f,   0.0f,  0.0f,  1.0f,
            -0.5f,  0.5f,  0.5f,  1.0f, 1.0f, 0.0f,   0.0f,  0.0f,  1.0f,
            -0.5f, -0.5f,  0.5f,  1.0f, 1.0f, 0.0f,   0.0f,  0.0f,  1.0f,
                                                     
            -0.5f,  0.5f,  0.5f,  0.0f, 0.0f, 1.0f,  -1.0f,  0.0f,  0.0f,
            -0.5f,  0.5f, -0.5f,  0.0f, 0.0f, 1.0f,  -1.0f,  0.0f,  0.0f,
            -0.5f, -0.5f, -0.5f,  0.0f, 0.0f, 1.0f,  -1.0f,  0.0f,  0.0f,
            -0.5f, -0.5f, -0.5f,  0.0f, 0.0f, 1.0f,  -1.0f,  0.0f,  0.0f,
            -0.5f, -0.5f,  0.5f,  0.0f, 0.0f, 1.0f,  -1.0f,  0.0f,  0.0f,
            -0.5f,  0.5f,  0.5f,  0.0f, 0.0f, 1.0f,  -1.0f,  0.0f,  0.0f,
                                                     
             0.5f,  0.5f,  0.5f,  0.0f, 1.0f, 0.0f,   1.0f,  0.0f,  0.0f,
             0.5f,  0.5f, -0.5f,  0.0f, 1.0f, 0.0f,   1.0f,  0.0f,  0.0f,
             0.5f, -0.5f, -0.5f,  0.0f, 1.0f, 0.0f,   1.0f,  0.0f,  0.0f,
             0.5f, -0.5f, -0.5f,  0.0f, 1.0f, 0.0f,   1.0f,  0.0f,  0.0f,
             0.5f, -0.5f,  0.5f,  0.0f, 1.0f, 0.0f,   1.0f,  0.0f,  0.0f,
             0.5f,  0.5f,  0.5f,  0.0f, 1.0f, 0.0f,   1.0f,  0.0f,  0.0f,
                                                     
            -0.5f, -0.5f, -0.5f,  1.0f, 1.0f, 1.0f,   0.0f, -1.0f,  0.0f,
             0.5f, -0.5f, -0.5f,  1.0f, 1.0f, 1.0f,   0.0f, -1.0f,  0.0f,
             0.5f, -0.5f,  0.5f,  1.0f, 1.0f, 1.0f,   0.0f, -1.0f,  0.0f,
             0.5f, -0.5f,  0.5f,  1.0f, 1.0f, 1.0f,   0.0f, -1.0f,  0.0f,
            -0.5f, -0.5f,  0.5f,  1.0f, 1.0f, 1.0f,   0.0f, -1.0f,  0.0f,
            -0.5f, -0.5f, -0.5f,  1.0f, 1.0f, 1.0f,   0.0f, -1.0f,  0.0f,
                                                     
            -0.5f,  0.5f, -0.5f,  0.0f, 1.0f, 1.0f,   0.0f,  1.0f,  0.0f,
             0.5f,  0.5f, -0.5f,  0.0f, 1.0f, 1.0f,   0.0f,  1.0f,  0.0f,
             0.5f,  0.5f,  0.5f,  0.0f, 1.0f, 1.0f,   0.0f,  1.0f,  0.0f,
             0.5f,  0.5f,  0.5f,  0.0f, 1.0f, 1.0f,   0.0f,  1.0f,  0.0f,
            -0.5f,  0.5f,  0.5f,  0.0f, 1.0f, 1.0f,   0.0f,  1.0f,  0.0f,
            -0.5f,  0.5f, -0.5f,  0.0f, 1.0f, 1.0f,   0.0f,  1.0f,  0.0f

        };
        // clang-format on

        std::vector<float> sphere_vertices = generateSphereVertices();

        int y = 11 % 9;
        std::cout << y;
        float max = 0.0f;
        for (int i = 0; i < sphere_vertices.size(); i++)
        {
            if (i % 9 < 3)
            {
                if (abs(sphere_vertices[i]) > max)
                {
                    max = abs(sphere_vertices[i]);
                }
            }

        }

        for (int i = 0; i < sphere_vertices.size(); i++)
        {
            if (i % 9 < 3)
            {
                sphere_vertices[i] /= max*2;
            }

        }

        //std::vector<float> vs = sphere_vertices;

        Program prog{ 1280, 720, "CSC305 A3" };
        Triangle cube{};

        cube.loadShaders();
        cube.loadDataToGPU(vertices);

        prog.run(cube);

        prog.freeGPUData();
        cube.freeGPUData();
    }
    catch (OpenGLError & err)
    {
        fmt::print("OpenGL Error:\n\t{}\n", err.what());
    }

    return 0;
}
