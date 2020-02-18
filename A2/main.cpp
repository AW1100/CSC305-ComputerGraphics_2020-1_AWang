#include "lab.hpp"

int main()
{
    using atlas::math::Point;
    using atlas::math::Ray;
    using atlas::math::Vector;

    std::shared_ptr<World> world{ std::make_shared<World>() };

    world->width = 600;
    world->height = 600;
    world->background = { 0, 0, 0 };
    world->sampler = std::make_shared<Regular>(4, 83);

    world->scene.push_back(
        std::make_shared<Sphere>(atlas::math::Point{ 50, 100, 0 }, 100.0f));
    world->scene[0]->setMaterial(
        std::make_shared<Matte>(0.50f, 0.05f, Colour{ 1, 0, 0 }));
    world->scene[0]->setColour({ 1, 0, 0 });

    world->scene.push_back(
        std::make_shared<Sphere>(atlas::math::Point{ 150, 150, 0 }, 50.0f));
    world->scene[1]->setMaterial(
        std::make_shared<Matte>(0.50f, 0.05f, Colour{ 0, 0, 1 }));
    world->scene[1]->setColour({ 0, 0, 1 });

    world->scene.push_back(
        std::make_shared<Triangle>(atlas::math::Point{ -200,200,50 }, atlas::math::Point{ -75,0,50 }, atlas::math::Point{ 0, 200, -200 }, atlas::math::Normal{ 0,0,1 }));
    world->scene[2]->setMaterial(
        std::make_shared<Matte>(0.50f, 0.05f, Colour{ 0, 1, 1 }));
    world->scene[2]->setColour({ 0, 1, 1 });

    world->scene.push_back(
        std::make_shared<Plane>(atlas::math::Point{ 0, 200, 0 }, atlas::math::Normal{ 0,-1,0 }));
    world->scene[3]->setMaterial(
        std::make_shared<Matte>(0.50f, 0.05f, Colour{ 0.75, 0.75, 0.75 }));
    world->scene[3]->setColour({ 0.75, 0.75, 0.75 });

    /*world.scene.push_back(
        std::make_shared<Sphere>(atlas::math::Point{ -128, 32, -700 }, 64.0f));
    world.scene[2]->setMaterial(
        std::make_shared<Matte>(0.50f, 0.05f, Colour{ 0, 1, 0 }));
    world.scene[2]->setColour({ 0, 1, 0 });*/

    world->ambient = std::make_shared<Ambient>();
    /*world->lights.push_back(
        std::make_shared<Directional>(Directional{ {0, 0, 1024} }));*/

    world->lights.push_back(
        std::make_shared<PointLight>(PointLight{ {200, -200, 100} }));

    world->ambient->setColour({ 1, 1, 1 });
    world->ambient->scaleRadiance(0.05f);

    world->lights[0]->setColour({ 1, 1, 1 });
    world->lights[0]->scaleRadiance(4.0f);

    /*world->lights[1]->setColour({ 1, 1, 1 });
    world->lights[1]->scaleRadiance(4.0f);*/

    Point samplePoint{}, pixelPoint{};

    Pinhole camera{};

    // change camera position here
    camera.setEye({ 0.0f, -200.0f, 1000.0f });

    camera.computeUVW();

    camera.renderScene(world);

    saveToFile("F:/csc305_spring2020_labs/labs/lab04_shading/raytrace.bmp", world->width, world->height, world->image);

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
