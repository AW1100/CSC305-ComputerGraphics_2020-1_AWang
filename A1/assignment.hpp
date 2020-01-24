#pragma once

#include <atlas/core/Float.hpp>
#include <atlas/math/Math.hpp>
#include <atlas/math/Ray.hpp>

#include <fmt/printf.h>
#include <stb_image.h>
#include <stb_image_write.h>

#include <vector>

#include <limits>

using Colour = atlas::math::Vector;

const std::size_t image_width{ 600 };
const std::size_t image_height{ 600 };
constexpr Colour background{ 0, 0, 0 };

void saveToBMP(std::string const& filename,
               std::size_t width,
               std::size_t height,
               std::vector<Colour> const& image);

// Your code here.

struct ShadeRec
{
    Colour colour;
    Colour cloest;
    float t;
};


class Sphere
{
public:
    constexpr Sphere() :
        center_{ 0,0,0 },
        radius_{ 0 },
        radius_squr_{ 0 },
        colour_{ 0,0,0 }
    {}

    constexpr Sphere(atlas::math::Point center, float radius, Colour colour) :
        center_{ center },
        radius_{ radius },
        radius_squr_{ radius * radius },
        colour_{ colour }
    {}

    bool hit(atlas::math::Ray<atlas::math::Vector> const& ray, ShadeRec& trace_data) const
    {
        auto o_c{ ray.o - center_ };

        auto a{ glm::dot(ray.d, ray.d) };
        auto b{ glm::dot(ray.d, o_c) * 2 };
        auto c{ glm::dot(o_c,o_c) - radius_squr_ };

        auto roots{ b * b - (4.0f * a * c) };

        if (roots >= 0.0f)
        {
            // hit
            trace_data.t = ((-b - std::sqrt(roots)) / (2.0f * a));


            trace_data.colour = colour_;

            return true;
        }
        return false;
    }

    Colour get_colour()
    {
        return colour_;
    }

    

private:
    atlas::math::Point center_;
    float radius_, radius_squr_;
    Colour colour_;

};

class Plane
{
public:
    constexpr Plane() :
        point_{ 0,0,0 },
        normal_{ 0,0,0 },
        colour_{ 0,0,0 }
    {}
    constexpr Plane(atlas::math::Point point, atlas::math::Normal normal, Colour colour) :
        point_{ point },
        normal_{ normal },
        colour_{ colour }
    {}

    bool hit(atlas::math::Ray<atlas::math::Vector> const& ray, ShadeRec& trace_data) const
    {
        //auto t = (point_ - ray.o) * normal_ / (ray.d * normal_);
        auto t = glm::dot((point_ - ray.o), normal_) / glm::dot(ray.d, normal_);

        if (t >= 0.0f)
        {
            trace_data.t = t;
            trace_data.colour = colour_;

            return true;
        }

        return false;
    }
    
private:
    atlas::math::Point point_;
    atlas::math::Normal normal_;
    Colour colour_;
};

//std::vector<Colour> sceneGenerator(std::vector<Colour> image, atlas::math::Ray<atlas::math::Vector> ray, Sphere object, ShadeRec& trace_data)
//{
//    for (std::size_t y{ 0 }; y < imageHeight; ++y)
//    {
//        for (std::size_t x{ 0 }; x < imageWidth; ++x)
//        {
//            // Compute origin of ray.
//            float originX = (x - 0.5f * (imageWidth - 1.0f));
//            float originY = (y - 0.5f * (imageHeight - 1.0f));
//
//            ray.o = { originX, originY, 100.0f };
//            
//            if (!object.hit(ray, trace_data))
//            {
//                trace_data.colour = background;
//            }
//
//            image[x + y * imageHeight] = trace_data.colour;
//        }
//    }
//
//    return image;
//};