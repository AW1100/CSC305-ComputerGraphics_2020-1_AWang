#include "assignment.hpp"
#include <iostream>
#include <fstream>
#include <string.h>
#include <string>
#include <sstream>

atlas::math::Point rotateYV(atlas::math::Point p, float angle);
atlas::math::Point rotateYN(atlas::math::Normal p, float angle);
atlas::math::Point scaleXYZ(atlas::math::Point p, float t);

int main()
{
    //std::shared_ptr<Mesh> mesh{ std::make_shared<Mesh>() };
    Mesh mesh{};
    //mesh.vertices = std::vector<atlas::math::Point>();
    std::ifstream myfile("../deer.obj");
    if (myfile.is_open())
    {
        std::string line;
        while (std::getline(myfile, line))
        {
            if (line[0] == 'v' && line[1] == ' ')
            {
                std::stringstream ls(line);
                std::string data;
                float vx, vy, vz;
                std::getline(ls, data, ' ');
                ls >> vx >> vy >> vz;
                atlas::math::Vector temp{ vx,-vy + 200,vz };
                mesh.vertices.push_back(temp);
            }

            if (line[0] == 'v' && line[1] == 'n')
            {
                std::stringstream ls(line);
                std::string data;
                float nx, ny, nz;
                std::getline(ls, data, ' ');
                ls >> nx >> ny >> nz;
                atlas::math::Normal temp{ nx,ny,nz };
                mesh.normals.push_back(temp);
            }

            if (line[0] == 'f' && line[1] == ' ')
            {
                std::stringstream ls(line);
                std::string data;
                std::string s1, s2, s3;
                std::getline(ls, data, ' ');
                ls >> s1 >> s2 >> s3;

                std::stringstream str1(s1);
                std::string sf1;
                int f1;
                std::getline(str1, sf1, '/');
                f1 = std::stoi(sf1);
                std::getline(str1, sf1, '/');
                int fn1;
                str1 >> fn1;


                std::stringstream str2(s2);
                std::string sf2;
                int f2;
                std::getline(str2, sf2, '/');
                f2 = std::stoi(sf2);
                std::getline(str2, sf2, '/');
                int fn2;
                str2 >> fn2;

                std::stringstream str3(s3);
                std::string sf3;
                int f3;
                std::getline(str3, sf3, '/');
                f3 = std::stoi(sf3);
                std::getline(str3, sf3, '/');
                int fn3;
                str2 >> fn3;

                std::vector vec{ f1,f2,f3 };
                mesh.vertex_faces.push_back(vec);
                std::vector vec2{ fn1,fn2,fn3 };
                mesh.vertex_norms.push_back(vec2);
            }
        }
    }
    myfile.close();

    mesh.num_vertices = (int)mesh.vertices.size();
    mesh.num_triangles = (int)mesh.vertex_faces.size();

    for (int i = 0; i < mesh.vertices.size(); i++)
    {
        mesh.vertices[i] = rotateYV(mesh.vertices[i], 45.0f);
        //mesh.vertices[i] = scaleXYZ(mesh.vertices[i], 100.0f);
    }

    for (int i = 0; i < mesh.normals.size(); i++)
    {
        mesh.normals[i] = rotateYN(mesh.normals[i], 45.0f);
    }


    std::shared_ptr<World> world{ std::make_shared<World>() };

    world->viewplane = std::make_shared<ViewPlane>(1000, 1000);
    world->background = { 0, 0, 0 };
    world->sampler = std::make_shared<Regular>(16, 83);

    world->objects.push_back(
        std::make_shared<Sphere>(atlas::math::Point{ 250, 100, 0 }, 100.0f));
    world->objects[0]->setMaterial(
        std::make_shared<Matte>(0.50f, 0.05f, Colour{ 1, 0, 0 }));
    world->objects[0]->setColour({ 1, 0, 0 });

    world->objects.push_back(
        std::make_shared<Sphere>(atlas::math::Point{ 380, 150, 0 }, 50.0f));
    world->objects[1]->setMaterial(
        std::make_shared<Matte>(0.50f, 0.05f, Colour{ 0, 0, 1 }));
    world->objects[1]->setColour({ 0, 0, 1 });

    world->objects.push_back(
        std::make_shared<Triangle>(atlas::math::Point{ -300,200,50 }, atlas::math::Point{ -175,0,50 }, atlas::math::Point{ -100, 200, -200 }, atlas::math::Normal{ 0,0,1 }));
    world->objects[2]->setMaterial(
        std::make_shared<Matte>(0.50f, 0.05f, Colour{ 0, 1, 1 }));
    world->objects[2]->setColour({ 0, 1, 1 });

    world->objects.push_back(
        std::make_shared<Box>(atlas::math::Point{ -500, 50, -200 }, atlas::math::Point{ -350, 200, -50 }));
    world->objects[3]->setMaterial(
        std::make_shared<Matte>(0.50f, 0.05f, Colour{ 1,1,0 }));
    world->objects[3]->setColour({ 1,1,0 });

    world->objects.push_back(
        std::make_shared<Plane>(atlas::math::Point{ 0, 200, 0 }, atlas::math::Normal{ 0,-1,0 }));
    world->objects[4]->setMaterial(
        std::make_shared<Matte>(0.50f, 0.05f, Colour{ 0.75, 0.75, 0.75 }));
    world->objects[4]->setColour({ 0.75, 0.75, 0.75 });

    world->objects.push_back(
        std::make_shared<Triangle>(atlas::math::Point{ -800,200,500 }, atlas::math::Point{ -800,-1000,-500 }, atlas::math::Point{ -800, 200, -500 }, atlas::math::Normal{ 1,0,0 }));
    world->objects[5]->setMaterial(
        std::make_shared<Matte>(0.50f, 0.05f, Colour{ 0.75, 0.75, 0.75 }));
    world->objects[5]->setColour({ 0.75, 0.75, 0.75 });

    world->objects.push_back(
        std::make_shared<Triangle>(atlas::math::Point{ 800,200,500 }, atlas::math::Point{ 800,-1000,-500 }, atlas::math::Point{ 800, 200, -500 }, atlas::math::Normal{ -1,0,0 }));
    world->objects[6]->setMaterial(
        std::make_shared<Matte>(0.50f, 0.05f, Colour{ 0.75, 0.75, 0.75 }));
    world->objects[6]->setColour({ 0.75, 0.75, 0.75 });

    world->objects.push_back(
        std::make_shared<Triangle>(atlas::math::Point{ -800,200,500 }, atlas::math::Point{ -800,-1000,-500 }, atlas::math::Point{ -800,-1000,500 }, atlas::math::Normal{ 1,0,0 }));
    world->objects[7]->setMaterial(
        std::make_shared<Matte>(0.50f, 0.05f, Colour{ 0.75, 0.75, 0.75 }));
    world->objects[7]->setColour({ 0.75, 0.75, 0.75 });

    world->objects.push_back(
        std::make_shared<Triangle>(atlas::math::Point{ 800,200,500 }, atlas::math::Point{ 800,-1000,-500 }, atlas::math::Point{ 800,-1000,500 }, atlas::math::Normal{ -1,0,0 }));
    world->objects[8]->setMaterial(
        std::make_shared<Matte>(0.50f, 0.05f, Colour{ 0.75, 0.75, 0.75 }));
    world->objects[8]->setColour({ 0.75, 0.75, 0.75 });

    world->objects.push_back(
        std::make_shared<Plane>(atlas::math::Point{ 0, 0, -500 }, atlas::math::Normal{ 0,0,1 }));
    world->objects[9]->setMaterial(
        std::make_shared<Matte>(0.50f, 0.05f, Colour{ 0.75, 0.75, 0.75 }));
    world->objects[9]->setColour({ 0.75, 0.75, 0.75 });

    world->objects.push_back(
        std::make_shared<Plane>(atlas::math::Point{ 0, -1000, 0 }, atlas::math::Normal{ 0,1,0 }));
    world->objects[10]->setMaterial(
        std::make_shared<Matte>(0.50f, 0.05f, Colour{ 0.75, 0.75, 0.75 }));
    world->objects[10]->setColour({ 0.75, 0.75, 0.75 });

    /*world->objects.push_back(
        std::make_shared<Sphere>(atlas::math::Point{ 0, 100, 0 }, 100.0f));
    world->objects[11]->setMaterial(
        std::make_shared<Matte>(0.50f, 0.05f, Colour{ 1, 0.5, 0.5 }));
    world->objects[11]->setColour({ 1, 0.5, 0.5 });*/

    for (int i = 11; i < mesh.num_triangles+11; i++)
    {
        world->objects.push_back(
            std::make_shared<MeshTriangle>(&mesh, mesh.vertex_faces[i-11][0]-1, mesh.vertex_faces[i-11][1]-1, mesh.vertex_faces[i-11][2]-1, mesh.normals[mesh.vertex_norms[i-11][0]-1]));
        world->objects[i]->setMaterial(
            std::make_shared<Matte>(0.50f, 0.05f, Colour{ 1, 0, 1 }));
        world->objects[i]->setColour({ 1, 0, 1 });
    }

    std::shared_ptr<Sampler> sampler2{ std::make_shared<Regular>(16, 83) };
    world->ambient = std::make_shared<AmbientOccluder>(sampler2);

    //world->ambient = std::make_shared<Ambient>();
    world->lights.push_back(
        std::make_shared<Directional>(Directional{ {0, 0, 1024} }));

    world->lights.push_back(
        std::make_shared<PointLight>(PointLight{ {500, -450, 200} }));

    /*world->lights.push_back(
        std::make_shared<AmbientOccluder>(AmbientOccluder{ sampler2 }));*/

    world->ambient->setColour({ 1, 1, 1 });
    world->ambient->scaleRadiance(1.0f);
    //world->ambient->setShadow(true);

    /*world->lights[0]->setColour({ 1, 1, 1 });
    world->lights[0]->scaleRadiance(6.0f);*/

    world->lights[0]->setColour({ 1, 1, 1 });
    world->lights[0]->scaleRadiance(2.0f);

    world->lights[1]->setColour({ 1, 1, 1 });
    world->lights[1]->scaleRadiance(5.0f);
    world->lights[1]->setShadow(true);

    /*world->lights[2]->setColour({ 1, 1, 1 });
    world->lights[2]->scaleRadiance(2.0f);*/
    //world->lights[2]->setShadow(true);

    Point3D samplePoint{}, pixelPoint{};

    Pinhole camera{};
    camera.setEye({ 0.0f, -200.0f, 1000.0f });

    camera.computeUVW();

    camera.renderScene(world);

    saveToFile("raytrace_Pinhole.bmp", world->viewplane->hres, world->viewplane->vres, world->image);


    return 0;
}

void saveToFile(std::string const& filename,
    std::size_t width,
    std::size_t height,
    std::vector<Colour> const& image)
{
    std::vector<unsigned char> data(image.size() * 3);

    for (std::size_t i{ 0 }, k{ 0 }; i < image.size(); ++i, k += 3)
    {
        Colour pixel = image[i];
        data[k + 0] = static_cast<unsigned char>(pixel.r * 255);
        data[k + 1] = static_cast<unsigned char>(pixel.g * 255);
        data[k + 2] = static_cast<unsigned char>(pixel.b * 255);
    }

    stbi_write_bmp(filename.c_str(),
        static_cast<int>(width),
        static_cast<int>(height),
        3,
        data.data());
}

atlas::math::Point rotateYV(atlas::math::Point p, float angle)
{

    angle = angle * 3.141592f / 180;
    float temp = p.z;
    p.z = p.z * cos(angle) - p.x * sin(angle);
    p.x = temp * sin(angle) + p.x * cos(angle);
    return p;
}

atlas::math::Point rotateYN(atlas::math::Normal p, float angle)
{

    angle = angle * 3.141592f / 180;
    float temp = p.z;
    p.z = p.z * cos(angle) - p.x * sin(angle);
    p.x = temp * sin(angle) + p.x * cos(angle);
    return p;
}

atlas::math::Point scaleXYZ(atlas::math::Point p, float t)
{
    p.x = t * p.x;
    p.y = t * p.y;
    p.z = t * p.z;
    return p;
}

