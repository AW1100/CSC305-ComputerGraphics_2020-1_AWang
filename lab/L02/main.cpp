#include "lab.hpp"

int main()
{
   // Your code here.
	constexpr std::size_t image_width { 512 };
	constexpr std::size_t image_height { 512 };
	constexpr Colour background{ 0, 0, 0 };
	constexpr Sphere object { {256,256,0},128,{1,0,0} };


	atlas::math::Ray<atlas::math::Vector> ray{ {0,0,0},{0,0,-1} };
	

	std::vector<Colour> image{ image_width * image_height };
	ShadeRec trace_data{};

	for (std::size_t y{ 0 }; y < image_height; y++)
	{
		for (std::size_t x{ 0 }; x < image_width; x++)
		{
			ray.o = { x + 0.5f, y + 0.5f, 0 };

			// did I hit
			if (!object.hit(ray, trace_data))
			{
				trace_data.colour = background;
			}

			image[x + y * image_height] = trace_data.colour;
		}
	}

	saveToFile("H:/csc305_spring2020_labs/labs/lab02_sphere/op.bmp", image_width, image_height, image);
	
    return 0;
}

void saveToFile(std::string const& filename,
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
