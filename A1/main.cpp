#include "assignment.hpp"

int main()
{
    // Your code here.
    atlas::math::Ray<atlas::math::Vector> ray{ {0,0,0},{0,0,-1} };
    Sphere spheres[num_sphere];
    spheres[0] = { {-50,0,0},60,{1,0,0} };
    spheres[1] = { {25,0,20},70,{0,1,0} };
    spheres[2] = { {-100,80,30},50,{0,0,1} };
    spheres[3] = { {0,90,0},50,{1,1,0} };
    spheres[4] = { {50,150,20},40,{1,0,1} };
    spheres[5] = { {100,50,0},90,{0,1,1} };
    Plane planes[num_plane];
    planes[0] = { { 200,0,0 }, { 1,0,1 }, { 0.75,0.5,0.25 } };
    planes[1] = { {-200,0,0},{0,1,1},{0.25,0.5,0.75} };
    Triangle triangles[num_triangle];
    triangles[0] = { {-150, -150, 0}, { 0,-250,0 }, { 20,-150,0 }, { 0,0.5,1 }, { 1,1,1 } };
    triangles[1] = { {-150, 0, 40}, { 0,250,30 }, { 50,150,20 }, { 0,0.5,1 }, { 0.5,1,1 } };

    std::vector<Colour> image(image_width * image_height);
    ShadeRec trace_data{};
    Sampler sampler{ 16,background,{0.0f,0.0f} };

    /*sceneGenerator(ray, trace_data, image, spheres, planes, triangles);
    saveToBMP("Finalscene.bmp", image_width, image_height, image);*/

    sceneGeneratorWithRegularMSAA(ray, trace_data, image, spheres, planes, triangles, sampler);
    saveToBMP("FinalsceneWithRegularMSAA.bmp", image_width, image_height, image);

    /*sceneGeneratorWithJitteredMSAA(ray, trace_data, image, spheres, planes, triangles, sampler);
    saveToBMP("FinalsceneWithJitteredMSAA.bmp", image_width, image_height, image);*/

    /*sceneGeneratorWithRandomMSAA(ray, trace_data, image, spheres, planes, triangles, sampler);
    saveToBMP("FinalsceneWithRandomMSAA.bmp", image_width, image_height, image);*/

}

/**
 * Saves a BMP image file based on the given array of pixels. All pixel values
 * have to be in the range [0, 1].
 *
 * @param filename The name of the file to save to.
 * @param width The width of the image.
 * @param height The height of the image.
 * @param image The array of pixels representing the image.
 */
void saveToBMP(std::string const& filename,
               std::size_t width,
               std::size_t height,
               std::vector<Colour> const& image)
{
    std::vector<unsigned char> data(image.size() * 3);

    for (std::size_t i{0}, k{0}; i < image.size(); ++i, k += 3)
    {
        Colour pixel = image[i];
        data[k + 0]  = static_cast<unsigned char>(pixel.r * 255);
        data[k + 1]  = static_cast<unsigned char>(pixel.g * 255);
        data[k + 2]  = static_cast<unsigned char>(pixel.b * 255);
    }

    stbi_write_bmp(filename.c_str(),
                   static_cast<int>(width),
                   static_cast<int>(height),
                   3,
                   data.data());
}
