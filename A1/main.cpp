#include "assignment.hpp"

int main()
{
    // Your code here.
    atlas::math::Ray<atlas::math::Vector> ray{ {0,0,0},{0,0,-1} };
    constexpr Sphere object{ {0,0,0},30,{1,0,0} };
    Sphere spheres[6];
    spheres[0] = { {-150,100,30},30,{1,0,0} };
    spheres[1] = { {-100,120,30},40,{0,1,0} };
    spheres[2] = { {-50,80,60},50,{0,0,1} };
    spheres[3] = { {0,90,190},60,{1,1,0} };
    spheres[4] = { {50,100,120},70,{1,0,1} };
    spheres[5] = { {100,50,150},80,{0,1,1} };
    Plane planes[2];
    planes[0] = { { 0,0,0 }, { 0,-1,-1 }, { 0.75,0.5,0.25 } };
    planes[1] = { {0,0,0},{1,0,1},{0.25,0.5,0.75} };

    std::vector<Colour> image(image_width * image_height);
    ShadeRec trace_data{};

    for (std::size_t y{ 0 }; y < image_height; y++)
    {
        for (std::size_t x{ 0 }; x < image_width; x++)
        {
            //ray.o = { x + 0.5f, y + 0.5f, 0 };
            float originX = (x - 0.5f * (image_width - 1.0f));
            float originY = (y - 0.5f * (image_height - 1.0f));

            ray.o = { originX, originY, 100.0f };

            float cloest_hit = std::numeric_limits<float>::infinity();
            // did I hit
            
            for (int i = 0; i < 6; i++)
            {
                if (!spheres[i].hit(ray, trace_data))
                {
                    if (cloest_hit < std::numeric_limits<float>::infinity())
                        continue;
                    trace_data.cloest = background;
                }
                else
                {
                    if (trace_data.t < cloest_hit)
                    {

                        cloest_hit = trace_data.t;
                        trace_data.cloest = trace_data.colour;
                    }
                    else
                    {
                        continue;
                    }
                    
                    //break;
                }
           
            }

            for (int i = 0; i < 2; i++)
            {
                if (planes[i].hit(ray, trace_data))
                {
                    if (trace_data.t < cloest_hit)
                    {
                        cloest_hit = trace_data.t;
                        trace_data.cloest = trace_data.colour;
                        
                    }

                }

            }


            /*if (trace_data.colour == background)
            {
                for (int i = 0; i < 2; i++)
                {
                    planes[i].hit(ray, trace_data);
                }
            }*/


            image[x + y * image_height] = trace_data.cloest;
        }
    }

    saveToBMP("Finalscene.bmp", image_width, image_height, image);
    
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
