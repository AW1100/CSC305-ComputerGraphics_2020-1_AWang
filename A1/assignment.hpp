#pragma once

#include <atlas/core/Float.hpp>
#include <atlas/math/Math.hpp>
#include <atlas/math/Ray.hpp>

#include <fmt/printf.h>
#include <stb_image.h>
#include <stb_image_write.h>

#include <vector>

using Colour = atlas::math::Vector;

void saveToBMP(std::string const& filename,
    std::size_t width,
    std::size_t height,
    std::vector<Colour> const& image);

// Your code here.
const int num_sphere{ 6 };
const int num_plane{ 2 };
const int num_triangle{ 2 };
const std::size_t image_width{ 600 };
const std::size_t image_height{ 600 };
constexpr Colour background{ 0, 0, 0 };

struct ShadeRec
{
    Colour colour;
    float t;
};

struct Sampler
{
    Sampler(int num_samples, Colour pixel_colour, atlas::math::Point2 pp) :
        num_samples_{ num_samples },
        pixel_colour_{ pixel_colour },
        n{ (int)sqrt((float)num_samples) },
        pp_{ pp }
    {}



    int num_samples_;
    Colour pixel_colour_;
    int n;
    atlas::math::Point2 pp_;
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

class Triangle
{
public:
    constexpr Triangle() :
        v0_{ 0,0,0 },
        v1_{ 0,0,0 },
        v2_{ 0,0,0 },
        normal_{ 0,0,0 },
        colour_{ background }


    {}
    constexpr Triangle(atlas::math::Point v0, atlas::math::Point v1, atlas::math::Point v2, atlas::math::Normal normal, Colour colour) :
        v0_{ v0 },
        v1_{ v1 },
        v2_{ v2 },
        normal_{ normal },
        colour_{ colour }


    {}

    // triangle hit function from "Ray Tracing From The Ground Up" textbook
    bool hit(atlas::math::Ray<atlas::math::Vector> const& ray, ShadeRec& trace_data) const
    {
        double a = v0_.x - v1_.x, b = v0_.x - v2_.x, c = ray.d.x, d = v0_.x - ray.o.x;
        double e = v0_.y - v1_.y, f = v0_.y - v2_.y, g = ray.d.y, h = v0_.y - ray.o.y;
        double i = v0_.z - v1_.z, j = v0_.z - v2_.z, k = ray.d.z, l = v0_.z - ray.o.z;

        double m = f * k - g * j, n = h * k - g * l, p = f * l - h * j;
        double q = g * i - e * k, s = e * j - f * i;

        double inv_denom = 1.0 / (a * m + b * q + c * s);

        double e1 = d * m - b * n - c * p;
        double beta = e1 * inv_denom;

        if (beta < 0.0)
            return (false);

        double r = e * l - h * i;
        double e2 = a * n + d * q + c * r;
        double gamma = e2 * inv_denom;

        if (gamma < 0.0)
            return (false);

        if (beta + gamma > 1.0)
            return (false);

        double e3 = a * p - b * r + d * s;
        double t = e3 * inv_denom;

        if (t < 0.0f)
            return (false);

        trace_data.t = (float)t;
        trace_data.colour = colour_;
        //sr.normal = normal;
        //sr.local_hit_point = ray.o + t * ray.d;

        return (true);
    }

private:
    atlas::math::Point v0_, v1_, v2_;
    atlas::math::Normal normal_;
    Colour colour_;
};

Colour objCasesDealing(atlas::math::Ray<atlas::math::Vector>& ray, ShadeRec& trace_data, Sphere spheres[], Plane planes[], Triangle triangles[])
{
    float cloest_t = -1.0f;
    Colour cloest_colour = background;
    for (int i = 0; i < num_sphere; i++)
    {
        if (!spheres[i].hit(ray, trace_data) && cloest_t < 0.0f)
        {
            trace_data.colour = background;

        }
        else if (!spheres[i].hit(ray, trace_data) && cloest_t >= 0.0f)
        {
            continue;
        }
        else if (spheres[i].hit(ray, trace_data) && cloest_t < 0.0f)
        {
            cloest_t = trace_data.t;
            cloest_colour = trace_data.colour;
        }
        else if (spheres[i].hit(ray, trace_data) && cloest_t >= 0.0f)
        {
            if (trace_data.t < cloest_t)
            {
                cloest_t = trace_data.t;
                cloest_colour = trace_data.colour;
            }
            else
            {
                continue;
            }
        }
    }

    for (int i = 0; i < num_plane; i++)
    {
        if (!planes[i].hit(ray, trace_data) && cloest_t < 0.0f)
        {
            trace_data.colour = background;

        }
        else if (!planes[i].hit(ray, trace_data) && cloest_t >= 0.0f)
        {
            continue;
        }
        else if (planes[i].hit(ray, trace_data) && cloest_t < 0.0f)
        {
            cloest_t = trace_data.t;
            cloest_colour = trace_data.colour;
        }
        else if (planes[i].hit(ray, trace_data) && cloest_t >= 0.0f)
        {
            if (trace_data.t < cloest_t)
            {
                cloest_t = trace_data.t;
                cloest_colour = trace_data.colour;
            }
            else
            {
                continue;
            }
        }
    }


    for (int i = 0; i < num_triangle; i++)
    {
        if (!triangles[i].hit(ray, trace_data) && cloest_t < 0.0f)
        {
            trace_data.colour = background;

        }
        else if (!triangles[i].hit(ray, trace_data) && cloest_t >= 0.0f)
        {
            continue;
        }
        else if (triangles[i].hit(ray, trace_data) && cloest_t < 0.0f)
        {
            cloest_t = trace_data.t;
            cloest_colour = trace_data.colour;
        }
        else if (triangles[i].hit(ray, trace_data) && cloest_t >= 0.0f)
        {
            if (trace_data.t < cloest_t)
            {
                cloest_t = trace_data.t;
                cloest_colour = trace_data.colour;
            }
            else
            {
                continue;
            }
        }
    }

    return cloest_colour;
}

void sceneGenerator(atlas::math::Ray<atlas::math::Vector>& ray, ShadeRec& trace_data, std::vector<Colour>& image, Sphere spheres[], Plane planes[], Triangle triangles[])
{
    


    for (std::size_t y{ 0 }; y < image_height; y++)
    {
        for (std::size_t x{ 0 }; x < image_width; x++)
        {
            float originX = (x - 0.5f * (image_width - 1.0f));
            float originY = (y - 0.5f * (image_height - 1.0f));

            ray.o = { originX, originY, 100.0f };           

            image[x + y * image_height] = objCasesDealing(ray, trace_data, spheres, planes, triangles);
        }
    }
};


void sceneGeneratorWithRegularMSAA(atlas::math::Ray<atlas::math::Vector>& ray, ShadeRec& trace_data, std::vector<Colour>& image, Sphere spheres[], Plane planes[], Triangle triangles[], Sampler& sampler)
{

    for (std::size_t y{ 0 }; y < image_height; y++)
    {
        for (std::size_t x{ 0 }; x < image_width; x++)
        {
            sampler.pixel_colour_ = background;

            for (int p = 0; p < sampler.n; p++)
            {
                for (int q = 0; q < sampler.n; q++)
                {
                    
                    sampler.pp_.x = (float)(1.0 * (x - 0.5 * image_width + (q + 0.5) / sampler.n));
                    sampler.pp_.y = (float)(1.0 * (y - 0.5 * image_height + (p + 0.5) / sampler.n));
                    ray.o = { sampler.pp_.x,sampler.pp_.y,100.0f };

                    sampler.pixel_colour_ += objCasesDealing(ray, trace_data, spheres, planes, triangles);
                }
            }
 
            sampler.pixel_colour_ /= sampler.num_samples_;
            image[x + y * image_height] = sampler.pixel_colour_;
        }
    }
};

void sceneGeneratorWithJitteredMSAA(atlas::math::Ray<atlas::math::Vector>& ray, ShadeRec& trace_data, std::vector<Colour>& image, Sphere spheres[], Plane planes[], Triangle triangles[], Sampler& sampler)
{

    for (std::size_t y{ 0 }; y < image_height; y++)
    {
        for (std::size_t x{ 0 }; x < image_width; x++)
        {
            sampler.pixel_colour_ = background;


            for (int p = 0; p < sampler.n; p++)
            {
                for (int q = 0; q < sampler.n; q++)
                {
                    
                    float r_x = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
                    float r_y = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
                    sampler.pp_.x = (float)(1.0 * (x - 0.5 * image_width + (q + r_x) / sampler.n));
                    sampler.pp_.y = (float)(1.0 * (y - 0.5 * image_height + (p + r_y) / sampler.n));
                    ray.o = { sampler.pp_.x,sampler.pp_.y,100.0f };

                    sampler.pixel_colour_ += objCasesDealing(ray, trace_data, spheres, planes, triangles);
                }
            }

            sampler.pixel_colour_ /= sampler.num_samples_;
            image[x + y * image_height] = sampler.pixel_colour_;
        }
    }
};

void sceneGeneratorWithRandomMSAA(atlas::math::Ray<atlas::math::Vector>& ray, ShadeRec& trace_data, std::vector<Colour>& image, Sphere spheres[], Plane planes[], Triangle triangles[], Sampler& sampler)
{

    for (std::size_t y{ 0 }; y < image_height; y++)
    {
        for (std::size_t x{ 0 }; x < image_width; x++)
        {
            sampler.pixel_colour_ = background;

            for (int j = 0; j < sampler.num_samples_; j++)
            {
                
                float r_x = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
                float r_y = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
                sampler.pp_.x = (float)(1.0 * (x - 0.5 * image_width + r_x));
                sampler.pp_.y = (float)(1.0 * (y - 0.5 * image_height + r_y));
                ray.o = { sampler.pp_.x,sampler.pp_.y,100.0f };

                sampler.pixel_colour_ += objCasesDealing(ray, trace_data, spheres, planes, triangles);
            }

            sampler.pixel_colour_ /= sampler.num_samples_;
            image[x + y * image_height] = sampler.pixel_colour_;
        }
    }
};

