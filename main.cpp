#include "assignment.hpp"

int main()
{
    // Your code here.
    const std::size_t imageWidth{ 600 };
    const std::size_t imageHeight{ 600 };

    Ray ray;
    ray.d = { 0, 0, -1 };
    Sphere s1;
    s1.init({ 125, 0, 0 }, 90.0f);
    Sphere s2;
    s2.init({ 75, 10, 0 }, 80.0f);
    Sphere s3;
    s3.init({ 25, 20, 0 }, 70.0f);
    Sphere s4;
    s4.init({ -25, 30, 0 }, 60.0f);
    Sphere s5;
    s5.init({ -75, 40, 0 }, 50.0f);
    Sphere s6;
    s6.init({ -125, 50, 0 }, 40.0f);
    //s.centre = { 0, 0, 0 };
    //s.radius = 120.0f;

    std::vector<Colour> image(imageWidth * imageHeight);
    //Colour pixel;
    /*
    for (std::size_t y{ 0 }; y < imageHeight; ++y)
    {
        for (std::size_t x{ 0 }; x < imageWidth; ++x)
        {
            // Compute origin of ray.
            float originX = (x - 0.5f * (imageWidth - 1.0f));
            float originY = (y - 0.5f * (imageHeight - 1.0f));

            ray.o = { originX, originY, 100.0f };

            pixel = intersectRayWithSphere(s1, ray, { 0,1,0 });


            image[x + y * imageHeight] = pixel;
        }
    }

    for (std::size_t y{ 0 }; y < imageHeight; ++y)
    {
        for (std::size_t x{ 0 }; x < imageWidth; ++x)
        {
            // Compute origin of ray.
            float originX = (x - 0.5f * (imageWidth - 1.0f));
            float originY = (y - 0.5f * (imageHeight - 1.0f));

            ray.o = { originX, originY, 100.0f };
            if (image[x + y * imageHeight] != black)
                continue;
            pixel = intersectRayWithSphere(s2, ray, { 1,0,0 });


            image[x + y * imageHeight] = pixel;
        }
    }
    */
    image = sceneGenerator(imageWidth, imageHeight, image, ray, s1, { 1,0,0 });
    image = sceneGenerator(imageWidth, imageHeight, image, ray, s2, { 0,1,0 });
    image = sceneGenerator(imageWidth, imageHeight, image, ray, s3, { 0,0,1 });
    image = sceneGenerator(imageWidth, imageHeight, image, ray, s4, { 1,1,0 });
    image = sceneGenerator(imageWidth, imageHeight, image, ray, s5, { 1,0,1 });
    image = sceneGenerator(imageWidth, imageHeight, image, ray, s6, { 0,1,1 });

    saveToBMP("sphere.bmp", imageWidth, imageHeight, image);
    
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
