#pragma once

#include <atlas/core/float.hpp>
#include <atlas/math/Math.hpp>
#include <atlas/math/Random.hpp>
#include <atlas/math/Ray.hpp>

#include <fmt/printf.h>
#include <stb_image.h>
#include <stb_image_write.h>

#include <limits>
#include <memory>
#include <vector>
#include <iostream>


using atlas::core::areEqual;

using Colour = atlas::math::Vector;
using Vector3D = atlas::math::Vector;
using Point3D = atlas::math::Point;
using Normal = atlas::math::Normal;
using Ray = atlas::math::Ray<atlas::math::Vector>;

void saveToFile(std::string const& filename,
    std::size_t width,
    std::size_t height,
    std::vector<Colour> const& image);

// Declarations
class BRDF;
class Camera;
class Material;
class Light;
class Shape;
class Sampler;

struct ViewPlane
{
    ViewPlane(int _hres, int _vres):
        hres{_hres},
        vres{_vres}
    {}

    int hres;
    int vres;
};

struct World
{
    std::shared_ptr<ViewPlane> viewplane;
    Colour background;
    std::shared_ptr<Sampler> sampler;
    std::vector<std::shared_ptr<Shape>> objects;
    std::vector<Colour> image;
    std::vector<std::shared_ptr<Light>> lights;
    std::shared_ptr<Light> ambient;
};

struct ShadeRec
{
    Colour colour;
    float t;
    Point3D hit_point;
    Point3D local_hit_point;
    Normal normal;
    Ray ray;
    std::shared_ptr<Material> material;
    std::shared_ptr<World> world;
};

class Sampler
{
public:
    Sampler(int numSamples, int numSets) :
        num_samples{ numSamples },
        num_sets{ numSets },
        count{ 0 },
        jump{ 0 }
    {
        samples.reserve(num_sets* num_samples);
        //hemisphere_samples.reserve(num_samples * num_sets);
        setupShuffledIndeces();
    }
    virtual ~Sampler() = default;

    int getNumSamples() const
    {
        return num_samples;
    }

    void setupShuffledIndeces()
    {
        shuffled_indices.reserve(num_samples * num_sets);
        std::vector<int> indices;

        std::random_device d;
        std::mt19937 generator(d());

        for (int j = 0; j < num_samples; ++j)
        {
            indices.push_back(j);
        }

        for (int p = 0; p < num_sets; ++p)
        {
            std::shuffle(indices.begin(), indices.end(), generator);

            for (int j = 0; j < num_samples; ++j)
            {
                shuffled_indices.push_back(indices[j]);
            }
        }
    }

    virtual void generateSamples() = 0;

    Point3D sampleUnitSquare()
    {
        if (count % num_samples == 0)
        {
            atlas::math::Random<int> engine;
            jump = (engine.getRandomMax() % num_sets) * num_samples;
        }

        return samples[jump + shuffled_indices[jump + count++ % num_samples]];
    }

    Point3D sampleHemisphere()
    {
        if (count % num_samples == 0) {	
            atlas::math::Random<int> engine;
            jump = (engine.getRandomMax() % num_sets) * num_samples;
        }
        return hemisphere_samples[jump + (count++) % num_samples];
    }

    void mapSamplesToHemisphere(const float exp) {
        int size = (int)samples.size();
        hemisphere_samples.reserve(num_samples * num_sets);

        for (int j = 0; j < size; j++) {
            float cos_phi = cos(2.0f * glm::pi<float>() * samples[j].x);
            float sin_phi = sin(2.0f * glm::pi<float>() * samples[j].x);
            float cos_theta = pow((1.0f - samples[j].y), 1.0f / (exp + 1.0f));
            float sin_theta = sqrt(1.0f - cos_theta * cos_theta);
            float pu = sin_theta * cos_phi;
            float pv = sin_theta * sin_phi;
            float pw = cos_theta;
            hemisphere_samples.push_back(Point3D(pu, pv, pw));
        }
    }

protected:
    std::vector<Point3D> samples;
    std::vector<Point3D> hemisphere_samples;
    std::vector<int> shuffled_indices;

    int num_samples;
    int num_sets;
    unsigned long count;
    int jump;
};

class Shape
{
public:
    Shape() : colour{ 0, 0, 0 }
    {}
    virtual ~Shape() = default;

    // if t computed is less than the t in sr, it and the colour should be
    // updated in sr
    virtual bool hit(Ray const& ray,
        ShadeRec& sr) const = 0;

    void setColour(Colour const& c)
    {
        colour = c;
    }

    Colour getColour() const
    {
        return colour;
    }

    void setMaterial(std::shared_ptr<Material> const& _material)
    {
        material = _material;
    }

    std::shared_ptr<Material> getMaterial() const
    {
        return material;
    }

    virtual bool shadowHit(Ray const& ray, float& tMin) const = 0;

    Point3D sample() {
        return Point3D(0.0f);
    }

    Normal getNormal() const {
        return Normal();
    }

protected:
    /*virtual bool intersectRay(Ray const& ray,
        float& tMin) const = 0;*/
    Colour colour;
    std::shared_ptr<Material> material;
};

class BRDF
{
public:
    virtual ~BRDF() = default;

    virtual Colour fn(ShadeRec const& sr,
        Vector3D const& reflected,
        Vector3D const& incoming) const = 0;
    virtual Colour rho(ShadeRec const& sr,
        Vector3D const& reflected) const = 0;
};

class Material
{
public:
    virtual ~Material() = default;

    virtual Colour shade(ShadeRec& sr) = 0;
};

class Light
{
public:
    Light():
        colour{ 1.0,1.0,1.0 },
        ls{ 1.0 },
        shadow{ true }
    {}

    virtual Vector3D getDirection(ShadeRec& sr) = 0;

    virtual Colour L([[maybe_unused]] ShadeRec& sr)
    {
        return ls * colour;
    }

    void scaleRadiance(float b)
    {
        ls = b;
    }

    void setColour(Colour const& c)
    {
        colour = c;
    }

    bool getShadow() const
    {
        return shadow;
    }

    void setShadow(const bool _shadow)
    {
        shadow = _shadow;
    }

    virtual bool inShadow(Ray& ray, ShadeRec& sr) const = 0;

protected:
    Colour colour;
    float ls;
    bool shadow;
};

class PointLight : public Light
{
public:
    PointLight() : Light{}
    {}

    PointLight(Point3D const& location) : Light{}
    {
        setLocation(location);
    }

    void setLocation(Vector3D const& _location)
    {
        location = _location;
    }

    Vector3D getDirection([[maybe_unused]] ShadeRec& sr)
    {
        return glm::normalize(location - sr.hit_point);
    }

    bool inShadow(Ray& ray, ShadeRec& sr) const
    {
        float t;
        int num_obj = (int)sr.world->objects.size();
        float d = sqrt((location.x - ray.o.x) * (location.x - ray.o.x) + (location.y - ray.o.y) * (location.y - ray.o.y) + (location.z - ray.o.z) * (location.z - ray.o.z));

        for (int i = 0; i < num_obj; i++)
        {
            if (sr.world->objects[i]->shadowHit(ray, t) && t < d)
                return true;
        }
        return false;

    }

private:
    Point3D location;
};

class Sphere : public Shape
{
public:
    Sphere(Point3D center, float radius) :
        center{ center }, radius{ radius }
    {}

    bool hit(Ray const& ray,
        ShadeRec& sr) const
    {
        Vector3D tmp = ray.o - center;
        float t{ std::numeric_limits<float>::max() };
        bool intersect{ intersectRay(ray, t) };

        // update ShadeRec info about new closest hit
        if (intersect && t < sr.t)
        {
            sr.normal = (tmp + t * ray.d) / radius;
            sr.ray = ray;
            sr.colour = colour;
            sr.t = t;
            sr.material = material;
        }

        return intersect;
    }
    bool shadowHit(Ray const& ray, float& tMin) const
    {
        float t;
        const auto tmp{ ray.o - center };
        const auto a{ glm::dot(ray.d, ray.d) };
        const auto b{ 2.0f * glm::dot(ray.d, tmp) };
        const auto c{ glm::dot(tmp, tmp) - radius * radius };
        const auto disc{ (b * b) - (4.0f * a * c) };
        const float kEpsilon{ 0.01f };
        if (disc < kEpsilon)
            return false;
        else {
            auto e = sqrt(disc);
            auto denom = 2.0 * a;
            t = (float)((-b - e) / denom);    // smaller root

            if (t > kEpsilon) {
                tMin = t;
                return true;
            }

            t = (float)((-b + e) / denom);    // larger root

            if (t > kEpsilon) {
                tMin = t;
                return true;
            }
        }
        return false;
    }

private:
    bool intersectRay(Ray const& ray,
        float& tMin) const
    {
        const auto tmp{ ray.o - center };
        const auto a{ glm::dot(ray.d, ray.d) };
        const auto b{ 2.0f * glm::dot(ray.d, tmp) };
        const auto c{ glm::dot(tmp, tmp) - radius * radius };
        const auto disc{ (b * b) - (4.0f * a * c) };

        if (atlas::core::geq(disc, 0.0f))
        {
            const float kEpsilon{ 0.01f };
            const float e{ std::sqrt(disc) };
            const float denom{ 2.0f * a };

            // Look at the negative root first
            float t = (-b - e) / denom;
            if (atlas::core::geq(t, kEpsilon))
            {
                tMin = t;
                return true;
            }

            // Now the positive root
            t = (-b + e);
            if (atlas::core::geq(t, kEpsilon))
            {
                tMin = t;
                return true;
            }
        }

        return false;
    }

    Point3D center;
    float radius;
};

class Box : public Shape
{
public:
    Box(Point3D _p0, Point3D _p1) :
        p0{ _p0 },
        p1{ _p1 },
        dims{ p1.x - p0.x,p1.y - p0.y,p1.z - p0.z }
    {}

    bool hit(Ray const& ray,
        ShadeRec& sr) const
    {
        float t{ std::numeric_limits<float>::max() };
        bool intersect{ intersectRay(ray, t, sr) };

        // update ShadeRec info about new closest hit
        if (intersect && t < sr.t)
        {
            //sr.normal = normal_;
            sr.ray = ray;
            sr.colour = colour;
            sr.t = t;
            sr.material = material;
        }

        return intersect;
    }

    bool shadowHit(Ray const& ray, float& tMin) const {
        Point3D o(ray.o);
        Point3D d(ray.d.x, ray.d.y, ray.d.z);
        Point3D t_min;
        Point3D t_max;

        float a = 1.0f / d.x;
        if (a >= 0.0f) {
            t_min.x = (p0.x - o.x) * a;
            t_max.x = (p1.x - o.x) * a;
        }
        else {
            t_min.x = (p1.x - o.x) * a;
            t_max.x = (p0.x - o.x) * a;
        }
        float b = 1.0f / d.y;
        if (b >= 0.0f) {
            t_min.y = (p0.y - o.y) * b;
            t_max.y = (p1.y - o.y) * b;
        }
        else {
            t_min.y = (p1.y - o.y) * b;
            t_max.y = (p0.y - o.y) * b;
        }
        float c = 1.0f / d.z;
        if (c >= 0.0f) {
            t_min.z = (p0.z - o.z) * c;
            t_max.z = (p1.z - o.z) * c;
        }
        else {
            t_min.z = (p1.z - o.z) * c;
            t_max.z = (p0.z - o.z) * c;
        }
        float t0, t1;
        int face_in, face_out;
        // finding largest
        if (t_min.x > t_min.y) {
            t0 = t_min.x;
            face_in = (a >= 0.0f) ? 0 : 3;
        }
        else {
            t0 = t_min.y;
            face_in = (b >= 0.0f) ? 1 : 4;
        }
        if (t_min.z > t0) {
            t0 = t_min.z;
            face_in = (c >= 0.0f) ? 2 : 5;
        }
        // find smallest
        if (t_max.x < t_max.y) {
            t1 = t_max.x;
            face_out = (a >= 0.0f) ? 3 : 0;
        }
        else {
            t1 = t_max.y;
            face_out = (b >= 0.0f) ? 4 : 1;
        }
        if (t_max.z < t1) {
            t1 = t_max.z;
            face_out = (c >= 0.0f) ? 5 : 2;
        }
        if (t0 < t1 && t1 > 0.01f) {
            if (t0 > 0.01f) {
                tMin = t0;
            }
            else {
                tMin = t1;
            }
            return true;
        }
        return false;
    }

    Normal get_normal(const int face_hit) const {
        switch (face_hit) {
        case 0: return (Normal(-1, 0, 0)); // -x face
        case 1: return (Normal(0, -1, 0)); // -y face
        case 2: return (Normal(0, 0, -1)); // -z face
        case 3: return (Normal(1, 0, 0));  // +x face
        case 4: return (Normal(0, 1, 0));  // +y face
        case 5: return (Normal(0, 0, 1));  // +z face
        }
        return (Normal(0, 0, 0));
    }

private:

    bool intersectRay(Ray const& ray,
        float& tMin, ShadeRec& sr) const
    {
        Point3D o(ray.o);
        Point3D d(ray.d.x, ray.d.y, ray.d.z);
        Point3D t_min;
        Point3D t_max;

        float a = 1.0f / d.x;
        if (a >= 0.0f) {
            t_min.x = (p0.x - o.x) * a;
            t_max.x = (p1.x - o.x) * a;
        }
        else {
            t_min.x = (p1.x - o.x) * a;
            t_max.x = (p0.x - o.x) * a;
        }
        float b = 1.0f / d.y;
        if (b >= 0.0f) {
            t_min.y = (p0.y - o.y) * b;
            t_max.y = (p1.y - o.y) * b;
        }
        else {
            t_min.y = (p1.y - o.y) * b;
            t_max.y = (p0.y - o.y) * b;
        }
        float c = 1.0f / d.z;
        if (c >= 0.0f) {
            t_min.z = (p0.z - o.z) * c;
            t_max.z = (p1.z - o.z) * c;
        }
        else {
            t_min.z = (p1.z - o.z) * c;
            t_max.z = (p0.z - o.z) * c;
        }
        float t0, t1;
        int face_in, face_out;
        // finding largest
        if (t_min.x > t_min.y) {
            t0 = t_min.x;
            face_in = (a >= 0.0f) ? 0 : 3;
        }
        else {
            t0 = t_min.y;
            face_in = (b >= 0.0f) ? 1 : 4;
        }
        if (t_min.z > t0) {
            t0 = t_min.z;
            face_in = (c >= 0.0f) ? 2 : 5;
        }
        // find smallest
        if (t_max.x < t_max.y) {
            t1 = t_max.x;
            face_out = (a >= 0.0f) ? 3 : 0;
        }
        else {
            t1 = t_max.y;
            face_out = (b >= 0.0f) ? 4 : 1;
        }
        if (t_max.z < t1) {
            t1 = t_max.z;
            face_out = (c >= 0.0f) ? 5 : 2;
        }
        if (t0 < t1 && t1 > 0.01f) {
            if (t0 > 0.01f) {
                tMin = t0;
                sr.normal = get_normal(face_in);
            }
            else {
                tMin = t1;
                sr.normal = get_normal(face_out);
            }
            //sr.hit_point = ray.o + sr.t * ray.d;
            return true;
        }
        else {
            return false;
        }
    }

    Point3D p0;
    Point3D p1;
    Point3D dims;
};

class Triangle : public Shape
{
public:

    Triangle(Point3D _v0, Point3D _v1, Point3D _v2, Normal _normal) :
        v0{ _v0 },
        v1{ _v1 },
        v2{ _v2 },
        normal{ _normal }

    {glm::normalize(normal); }

    // triangle hit function from "Ray Tracing From The Ground Up" textbook
    bool hit(Ray const& ray, ShadeRec& sr) const
    {
        float t{ std::numeric_limits<float>::max() };
        bool intersect{ intersectRay(ray, t) };

        // update ShadeRec info about new closest hit
        if (intersect && t < sr.t)
        {
            sr.normal = normal;
            sr.ray = ray;
            sr.colour = colour;
            sr.t = t;
            sr.material = material;
        }

        return intersect;
    }

    bool shadowHit(Ray const& ray, float& tMin) const
    {
        float a = v0.x - v1.x, b = v0.x - v2.x, c = ray.d.x, d = v0.x - ray.o.x;
        float e = v0.y - v1.y, f = v0.y - v2.y, g = ray.d.y, h = v0.y - ray.o.y;
        float i = v0.z - v1.z, j = v0.z - v2.z, k = ray.d.z, l = v0.z - ray.o.z;

        float m = f * k - g * j, n = h * k - g * l, p = f * l - h * j;
        float q = g * i - e * k, s = e * j - f * i;

        float inv_denom = 1.0f / (a * m + b * q + c * s);

        float e1 = d * m - b * n - c * p;
        float beta = e1 * inv_denom;
        const float kEpsilon{ 0.01f };
        if (beta < kEpsilon)
            return (false);

        float r = e * l - h * i;
        float e2 = a * n + d * q + c * r;
        float gamma = e2 * inv_denom;

        if (gamma < kEpsilon)
            return (false);

        if (beta + gamma > 1.0)
            return (false);

        float e3 = a * p - b * r + d * s;
        float t = e3 * inv_denom;

        if (t < kEpsilon)
            return (false);

        tMin = t;
        return (true);
    }

private:
    bool intersectRay(Ray const& ray,
        float& tMin) const
    {
        float a = v0.x - v1.x, b = v0.x - v2.x, c = ray.d.x, d = v0.x - ray.o.x;
        float e = v0.y - v1.y, f = v0.y - v2.y, g = ray.d.y, h = v0.y - ray.o.y;
        float i = v0.z - v1.z, j = v0.z - v2.z, k = ray.d.z, l = v0.z - ray.o.z;

        float m = f * k - g * j, n = h * k - g * l, p = f * l - h * j;
        float q = g * i - e * k, s = e * j - f * i;

        float inv_denom = 1.0f / (a * m + b * q + c * s);

        float e1 = d * m - b * n - c * p;
        float beta = e1 * inv_denom;

        if (beta < 0.0f)
            return (false);

        float r = e * l - h * i;
        float e2 = a * n + d * q + c * r;
        float gamma = e2 * inv_denom;

        if (gamma < 0.0f)
            return (false);

        if (beta + gamma > 1.0f)
            return (false);

        float e3 = a * p - b * r + d * s;
        float t = e3 * inv_denom;

        if (t < 0.0f)
            return (false);

        tMin = t;
        return (true);
    }

    Point3D v0, v1, v2;
    Normal normal;
};

class Plane : public Shape
{
public:

    Plane(Point3D _point, Normal _normal) :
        point{ _point },
        normal{ _normal }
    {glm::normalize(normal); }

    bool hit(Ray const& ray, ShadeRec& sr) const
    {
        float t{ std::numeric_limits<float>::max() };
        bool intersect{ intersectRay(ray, t) };

        // update ShadeRec info about new closest hit
        if (intersect && t < sr.t)
        {
            sr.normal = normal;
            sr.ray = ray;
            sr.colour = colour;
            sr.t = t;
            sr.material = material;
        }

        return intersect;
    }

    bool shadowHit(Ray const& ray, float& tMin) const
    {
        float t = glm::dot(point - ray.o, normal) / glm::dot(ray.d, normal);
        const float kEpsilon{ 0.01f };
        if (t > kEpsilon)
        {
            tMin = t;
            return true;
        }
        else
            return false;
    }


private:
    bool intersectRay(Ray const& ray,
        float& tMin) const
    {
        auto t = glm::dot((point - ray.o), normal) / glm::dot(ray.d, normal);

        const float kEpsilon{ 0.01f };
        if (atlas::core::geq(t, kEpsilon))
        {
            tMin = t;

            return true;
        }

        return false;
    }

    Point3D point;
    Normal normal;
};

class Regular : public Sampler
{
public:
    Regular(int num_samples, int num_sets) : Sampler{ num_samples, num_sets }
    {
        generateSamples();
    }

    void generateSamples()
    {
        int n = static_cast<int>(glm::sqrt(static_cast<float>(num_samples)));

        for (int j = 0; j < num_sets; ++j)
        {
            for (int p = 0; p < n; ++p)
            {
                for (int q = 0; q < n; ++q)
                {
                    samples.push_back(
                        Point3D{ (q + 0.5f) / n, (p + 0.5f) / n, 0.0f });
                }
            }
        }
    }
};

class Lambertian : public BRDF
{
public:
    Lambertian() : 
        cd{},
        kd{}
    {}
    Lambertian(float diffuseReflection, Colour diffuseColor) :
        kd{ diffuseReflection },
        cd{ diffuseColor }
    {}

    Colour
        fn([[maybe_unused]] ShadeRec const& sr,
            [[maybe_unused]] Vector3D const& reflected,
            [[maybe_unused]] Vector3D const& incoming) const
    {
        return kd * cd * glm::one_over_pi<float>();
    }

    Colour
        rho([[maybe_unused]] ShadeRec const& sr,
            [[maybe_unused]] Vector3D const& reflected) const
    {
        return cd * kd;
    }

    void setDiffuseReflection(float _kd)
    {
        kd = _kd;
    }

    void setDiffuseColour(Colour const& colour)
    {
        cd = colour;
    }

private:
    float kd;
    Colour cd;
};

class Matte : public Material
{
public:
    Matte() :
        Material{},
        ambient_brdf{ std::make_shared<Lambertian>() },
        diffuse_brdf{ std::make_shared<Lambertian>() }
    {}

    Matte(float kd, float ka, Colour colour) : Matte{}
    {
        setDiffuseReflection(kd);
        setAmbientReflection(ka);
        setDiffuseColour(colour);
    }

    void setDiffuseReflection(float k)
    {
        diffuse_brdf->setDiffuseReflection(k);
    }

    void setAmbientReflection(float k)
    {
        ambient_brdf->setDiffuseReflection(k);
    }

    void setDiffuseColour(Colour colour)
    {
        diffuse_brdf->setDiffuseColour(colour);
        ambient_brdf->setDiffuseColour(colour);
    }

    Colour shade(ShadeRec& sr)
    {
        Vector3D wo = -sr.ray.o;
        Colour L = ambient_brdf->rho(sr, wo) * sr.world->ambient->L(sr);
        size_t numLights = sr.world->lights.size();

        for (size_t i{ 0 }; i < numLights; ++i)
        {
            Vector3D wi = sr.world->lights[i]->getDirection(sr);
            float nDotWi = glm::dot(sr.normal, wi);

            if (nDotWi > 0.0f)
            {
                bool in_shadow = false;
                if (sr.world->lights[i]->getShadow())
                {
                    Ray shadowRay(sr.hit_point, wi);
                    in_shadow = sr.world->lights[i]->inShadow(shadowRay, sr);
                }

                if (!in_shadow)
                {
                    L += diffuse_brdf->fn(sr, wo, wi) * sr.world->lights[i]->L(sr) * nDotWi;
                }

            }
        }

        return L;
    }



private:
    std::shared_ptr<Lambertian> ambient_brdf;
    std::shared_ptr<Lambertian> diffuse_brdf;

};

class Directional : public Light
{
public:
    Directional() : Light{}
    {}

    Directional(Vector3D const& d) : Light{}
    {
        setDirection(d);
    }

    void setDirection(Vector3D const& d)
    {
        direction = glm::normalize(d);
    }

    Vector3D getDirection([[maybe_unused]] ShadeRec& sr)
    {
        return direction;
    }

    bool inShadow([[maybe_unused]] Ray& ray, [[maybe_unused]] ShadeRec& sr) const
    {
        return false;
    }

private:
    Vector3D direction;
};

class Ambient : public Light
{
public:
    Ambient() : Light{}
    {}

    Vector3D getDirection([[maybe_unused]] ShadeRec& sr)
    {
        return Vector3D{ 0.0f };
    }

    bool inShadow([[maybe_unused]] Ray& ray, [[maybe_unused]] ShadeRec& sr) const
    {
        return false;
    }

    Colour L([[maybe_unused]] ShadeRec& sr)
    {
        return ls * colour;
    }

};

class Camera
{
public:
    Camera() :
        eye{ 0.0f, 0.0f, 500.0f },
        lookAt{ 0.0f },
        up{ 0.0f, 1.0f, 0.0f },
        u{ 1.0f, 0.0f, 0.0f },
        v{ 0.0f, 1.0f, 0.0f },
        w{ 0.0f, 0.0f, 1.0f }
    {}

    virtual ~Camera() = default;

    virtual void renderScene(std::shared_ptr<World>& world) const = 0;

    void setEye(Point3D const& _eye)
    {
        eye = _eye;
    }

    void setLookAt(Point3D const& _lookAt)
    {
        lookAt = _lookAt;
    }

    void setUpVector(Vector3D const& _up)
    {
        up = _up;
    }

    void computeUVW()
    {
        w = glm::normalize(eye - lookAt);
        u = glm::normalize(glm::cross(up, w));
        v = glm::cross(w, u);

        if (areEqual(eye.x, lookAt.x) && areEqual(eye.z, lookAt.z) &&
            eye.y > lookAt.y)
        {
            u = { 0.0f, 0.0f, 1.0f };
            v = { 1.0f, 0.0f, 0.0f };
            w = { 0.0f, 1.0f, 0.0f };
        }

        if (areEqual(eye.x, lookAt.x) && areEqual(eye.z, lookAt.z) &&
            eye.y < lookAt.y)
        {
            u = { 1.0f, 0.0f, 0.0f };
            v = { 0.0f, 0.0f, 1.0f };
            w = { 0.0f, -1.0f, 0.0f };
        }
    }

protected:
    Point3D eye;
    Point3D lookAt;
    Point3D up;
    Vector3D u, v, w;
};

class Pinhole : public Camera
{
public:
    Pinhole() : Camera{}, distance{ 500.0f }, zoom{ 1.0f }
    {}

    void setDistance(float _distance)
    {
        distance = _distance;
    }

    void setZoom(float _zoom)
    {
        zoom = _zoom;
    }

    Vector3D rayDirection(Point3D const& p) const
    {
        const auto dir = p.x * u + p.y * v - distance * w;
        return glm::normalize(dir);
    }

    void renderScene(std::shared_ptr<World>& world) const
    {
        Point3D samplePoint{}, pixelPoint{};
        Ray ray{};

        ray.o = eye;
        float avg{ 1.0f / world->sampler->getNumSamples() };

        for (int r{ 0 }; r < world->viewplane->vres; ++r)
        {
            for (int c{ 0 }; c < world->viewplane->hres; ++c)
            {
                Colour pixelAverage{ 0, 0, 0 };

                for (int j = 0; j < world->sampler->getNumSamples(); ++j)
                {
                    ShadeRec trace_data{};
                    trace_data.world = world;
                    trace_data.t = std::numeric_limits<float>::max();
                    samplePoint = world->sampler->sampleUnitSquare();
                    pixelPoint.x = c - 0.5f * world->viewplane->hres + samplePoint.x;
                    pixelPoint.y = r - 0.5f * world->viewplane->vres + samplePoint.y;
                    ray.d = rayDirection(pixelPoint);

                    bool hit{};
                    for (auto obj : world->objects)
                    {
                        hit |= obj->hit(ray, trace_data);
                    }

                    if (hit)
                    {
                        trace_data.hit_point = ray.o + trace_data.t * ray.d;
                        pixelAverage += trace_data.material->shade(trace_data);
                    }
                }

                pixelAverage.r *= avg;
                pixelAverage.g *= avg;
                pixelAverage.b *= avg;
                // out-of-gamut handling Max-to-one
                float MAX_COLOUR = std::max(std::max(pixelAverage.r, pixelAverage.g), pixelAverage.b);
                if (MAX_COLOUR > 1.0f)
                {
                    pixelAverage.r /= MAX_COLOUR;
                    pixelAverage.g /= MAX_COLOUR;
                    pixelAverage.b /= MAX_COLOUR;
                }

                world->image.push_back({ pixelAverage.r,
                                       pixelAverage.g,
                                       pixelAverage.b });
            }
        }
    }

private:
    float distance;
    float zoom;
};

class Mesh
{
public:
    std::vector<Point3D> vertices;
    std::vector<Normal> normals;
    //std::vector<float> u;
    //std::vector<float> v;
    std::vector<std::vector<int>> vertex_faces;
    std::vector<std::vector<int>> vertex_norms;
    int num_vertices;
    int num_triangles;

};

class MeshTriangle : public Shape
{
public:
    Mesh* mesh_ptr;
    int index0, index1, index2;
    Normal normal;
    //float area;

    MeshTriangle(Mesh* _mesh_ptr, const int i0, const int i1, const int i2, Normal n) :Shape(),
        mesh_ptr{ _mesh_ptr },
        index0{ i0 }, index1{ i1 }, index2{ i2 },
        normal{ n }
    {}

    void computeNormal(const bool reverse_normal)
    {
        normal = glm::cross((mesh_ptr->vertices[index1] - mesh_ptr->vertices[index0]), (mesh_ptr->vertices[index2] - mesh_ptr->vertices[index0]));

        glm::normalize(normal);
        if (reverse_normal)
            normal = -normal;

    }

    bool hit(Ray const& ray, ShadeRec& sr) const
    {
        float t{ std::numeric_limits<float>::max() };
        bool intersect{ intersectRay(ray, t) };

        // update ShadeRec info about new closest hit
        if (intersect && t < sr.t)
        {
            sr.normal = normal;
            sr.ray = ray;
            sr.colour = colour;
            sr.t = t;
            sr.material = material;
        }

        return intersect;
    }

    bool shadowHit(Ray const& ray, float& tMin) const
    {
        float a = mesh_ptr->vertices[index0].x - mesh_ptr->vertices[index1].x, b = mesh_ptr->vertices[index0].x - mesh_ptr->vertices[index2].x, c = ray.d.x, d = mesh_ptr->vertices[index0].x - ray.o.x;
        float e = mesh_ptr->vertices[index0].y - mesh_ptr->vertices[index1].y, f = mesh_ptr->vertices[index0].y - mesh_ptr->vertices[index2].y, g = ray.d.y, h = mesh_ptr->vertices[index0].y - ray.o.y;
        float i = mesh_ptr->vertices[index0].z - mesh_ptr->vertices[index1].z, j = mesh_ptr->vertices[index0].z - mesh_ptr->vertices[index2].z, k = ray.d.z, l = mesh_ptr->vertices[index0].z - ray.o.z;

        float m = f * k - g * j, n = h * k - g * l, p = f * l - h * j;
        float q = g * i - e * k, s = e * j - f * i;

        float inv_denom = 1.0f / (a * m + b * q + c * s);

        float e1 = d * m - b * n - c * p;
        float beta = e1 * inv_denom;
        const float kEpsilon{ 0.01f };
        if (beta < kEpsilon)
            return (false);

        float r = e * l - h * i;
        float e2 = a * n + d * q + c * r;
        float gamma = e2 * inv_denom;

        if (gamma < kEpsilon)
            return (false);

        if (beta + gamma > 1.0)
            return (false);

        float e3 = a * p - b * r + d * s;
        float t = e3 * inv_denom;

        if (t < kEpsilon)
            return (false);

        tMin = t;
        return (true);
    }

private:

    bool intersectRay(Ray const& ray,
        float& tMin) const
    {
        float a = mesh_ptr->vertices[index0].x - mesh_ptr->vertices[index1].x, b = mesh_ptr->vertices[index0].x - mesh_ptr->vertices[index2].x, c = ray.d.x, d = mesh_ptr->vertices[index0].x - ray.o.x;
        float e = mesh_ptr->vertices[index0].y - mesh_ptr->vertices[index1].y, f = mesh_ptr->vertices[index0].y - mesh_ptr->vertices[index2].y, g = ray.d.y, h = mesh_ptr->vertices[index0].y - ray.o.y;
        float i = mesh_ptr->vertices[index0].z - mesh_ptr->vertices[index1].z, j = mesh_ptr->vertices[index0].z - mesh_ptr->vertices[index2].z, k = ray.d.z, l = mesh_ptr->vertices[index0].z - ray.o.z;

        float m = f * k - g * j, n = h * k - g * l, p = f * l - h * j;
        float q = g * i - e * k, s = e * j - f * i;

        float inv_denom = 1.0f / (a * m + b * q + c * s);

        float e1 = d * m - b * n - c * p;
        float beta = e1 * inv_denom;

        if (beta < 0.0f)
            return (false);

        float r = e * l - h * i;
        float e2 = a * n + d * q + c * r;
        float gamma = e2 * inv_denom;

        if (gamma < 0.0f)
            return (false);

        if (beta + gamma > 1.0f)
            return (false);

        float e3 = a * p - b * r + d * s;
        float t = e3 * inv_denom;

        if (t < 0.0f)
            return (false);

        tMin = t;
        return (true);
    }



};

class AmbientOccluder : public Light {
public:
    AmbientOccluder(std::shared_ptr<Sampler> sampler) :Light(),
        min_amount{0.25f},
        sampler_ptr{sampler}
    {sampler_ptr->mapSamplesToHemisphere(1.0f); }
    
    Vector3D getDirection([[maybe_unused]] ShadeRec& sr) {
        Point3D sp = sampler_ptr->sampleHemisphere();
        return (sp.x * u + sp.y * v + sp.z * w);
    }
    
    bool inShadow(Ray& ray, ShadeRec& sr) const
    {
        float t;
        int num_objects = (int)sr.world->objects.size();

        for (int j = 0; j < num_objects; j++)
            if (sr.world->objects[j]->shadowHit(ray, t))
                return true;

        return false;
    }

    Colour L(ShadeRec& sr) {
        w = glm::normalize(sr.normal);
        v = glm::normalize(glm::cross(w, Vector3D(0.0072f, 1.0f, 0.0034f))); // jitter the up vector in case normal is vertical
        //v = glm::normalize(glm::cross(w, Vector3D(0.0f, -1.0f, 0.0f)));
        u = glm::cross(v, w);

        Ray shadow_ray;
        shadow_ray.o = sr.hit_point;
        shadow_ray.d = getDirection(sr);

        if (inShadow(shadow_ray, sr))
        {
            return (min_amount * ls * colour);
        }
        else 
        {
            return (ls * colour);
        }
    }

private:
    float 		min_amount;
    Vector3D u, v, w;
    std::shared_ptr<Sampler> sampler_ptr;
};



