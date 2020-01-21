#pragma once

#include <atlas/math/Math.hpp>
#include <atlas/math/Ray.hpp>

#include <fmt/printf.h>
#include <stb_image.h>
#include <stb_image_write.h>

#include <vector>

using namespace atlas;
using Colour = math::Vector;

void saveToBMP(std::string const& filename,
               std::size_t width,
               std::size_t height,
               std::vector<Colour> const& image);

// Your code here.
static Colour background { 0, 0, 0 };

struct Ray
{
    atlas::math::Point o;
    atlas::math::Vector d;

    atlas::math::Point eval(float t)
    {
        return o + t * d;
    }
};

struct Sphere
{
    atlas::math::Point centre;
    float radius;

    void init(atlas::math::Point c, float r)
    {
        centre = c;
        radius = r;
    }
};

Colour intersectRayWithSphere(Sphere s, Ray r, Colour col)
{
    using atlas::math::Point;
    using atlas::math::Vector;

    Vector temp = r.o - s.centre;
    float a = glm::dot(r.d, r.d);
    float b = 2.0f * glm::dot(temp, r.d);
    float c = glm::dot(temp, temp) - (s.radius * s.radius);
    float disc = b * b - 4.0f * a * c;

    if (disc < 0.0f)
    {
        return background;
    }

    return col;
};

std::vector<Colour> sceneGenerator(std::size_t imageWidth, std::size_t imageHeight, std::vector<Colour> image, Ray ray, Sphere s, Colour col)
{
    Colour pixel;
    for (std::size_t y{ 0 }; y < imageHeight; ++y)
    {
        for (std::size_t x{ 0 }; x < imageWidth; ++x)
        {
            // Compute origin of ray.
            float originX = (x - 0.5f * (imageWidth - 1.0f));
            float originY = (y - 0.5f * (imageHeight - 1.0f));

            ray.o = { originX, originY, 100.0f };
            
            pixel = intersectRayWithSphere(s, ray, col);

            if (image[x + y * imageHeight] != background && pixel == background)
            {
                continue;
            }

            image[x + y * imageHeight] = pixel;
        }
    }

    return image;
};