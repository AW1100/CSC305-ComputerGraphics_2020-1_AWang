#pragma once

#include <atlas/core/Float.hpp>
#include <atlas/math/Math.hpp>
#include <atlas/math/Ray.hpp>

#include <fmt/printf.h>
#include <stb_image.h>
#include <stb_image_write.h>

#include <vector>

using Colour = atlas::math::Vector;

void saveToFile(std::string const& filename,
                std::size_t width,
                std::size_t height,
                std::vector<Colour> const& image);

// Your code here.

struct ShadeRec
{
	Colour colour;
	float t;
};

class Sphere
{
public:
	constexpr Sphere(atlas::math::Point center, float radius, Colour colour) :
		center_{ center },
		radius_{ radius },
		radius_squr_{radius * radius},
		colour_{colour}
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

private:
	atlas::math::Point center_;
		float radius_, radius_squr_;
	Colour colour_;
	
};

