#pragma once

#include <atlas/core/Float.hpp>
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

struct World
{
    std::size_t width, height;
    Colour background;
    std::shared_ptr<Sampler> sampler;
    std::vector<std::shared_ptr<Shape>> scene;
    std::vector<Colour> image;
    std::vector<std::shared_ptr<Light>> lights;
    std::shared_ptr<Light> ambient;
};

struct ShadeRec
{
    Colour color;
    float t;
    atlas::math::Point hitPoint;
    atlas::math::Normal normal;
    atlas::math::Ray<atlas::math::Vector> ray;
    std::shared_ptr<Material> material;
    std::shared_ptr<World> world;
};

// Abstract classes defining the interfaces for concrete entities

class Sampler
{
public:
    Sampler(int numSamples, int numSets) :
        mNumSamples{ numSamples }, mNumSets{ numSets }, mCount{ 0 }, mJump{ 0 }
    {
        mSamples.reserve(mNumSets* mNumSamples);
        //disk_samples.reserve(100);
        //disk_samples.reserve(mNumSets* mNumSamples);
        setupShuffledIndeces();
    }
    virtual ~Sampler() = default;

    int getNumSamples() const
    {
        return mNumSamples;
    }

    void setupShuffledIndeces()
    {
        mShuffledIndeces.reserve(mNumSamples * mNumSets);
        std::vector<int> indices;

        std::random_device d;
        std::mt19937 generator(d());

        for (int j = 0; j < mNumSamples; ++j)
        {
            indices.push_back(j);
        }

        for (int p = 0; p < mNumSets; ++p)
        {
            std::shuffle(indices.begin(), indices.end(), generator);

            for (int j = 0; j < mNumSamples; ++j)
            {
                mShuffledIndeces.push_back(indices[j]);
            }
        }
    }

    virtual void generateSamples() = 0;

    atlas::math::Point sampleUnitSquare()
    {
        if (mCount % mNumSamples == 0)
        {
            atlas::math::Random<int> engine;
            mJump = (engine.getRandomMax() % mNumSets) * mNumSamples;
        }

        return mSamples[mJump + mShuffledIndeces[mJump + mCount++ % mNumSamples]];
    }

    void map_samples_to_unit_disk()
    {
        //size_t size = mSamples.size();
        float r, phi;
        atlas::math::Point sp{};
        //disk_samples.reserve(4);
        for (size_t j = 0; j < 4; j++) {
            sp.x = (float)(2.0 * mSamples[j].x - 1.0);
            sp.y = (float)(2.0 * mSamples[j].y - 1.0);

            if (sp.x > -sp.y) {
                if (sp.x > sp.y) {
                    r = sp.x;
                    phi = sp.y / sp.x;
                }
                else {
                    r = sp.y;
                    phi = 2.0f - sp.x / sp.y;
                }
            }
            else {
                if (sp.x < sp.y) {
                    r = -sp.x;
                    phi = 4.0f + sp.y / sp.x;
                }
                else {
                    r = -sp.y;
                    if (sp.y != 0.0f) {
                        phi = 6.0f - sp.x / sp.y;
                    }
                    else {
                        phi = 0.0f;
                    }
                }
            }
            phi *= 3.14159f / 4.0f;

            disk_samples[j].x = r * cos(phi);
            disk_samples[j].y = r * sin(phi);
        }
        mSamples.erase(mSamples.begin(), mSamples.end());
    }

    atlas::math::Point sampleUnitDisk()
    {
        if (mCount % mNumSamples == 0)
        {
            atlas::math::Random<int> engine;
            mJump = (engine.getRandomMax() % mNumSets) * mNumSamples;
        }

        return disk_samples[mJump + mCount++ % mNumSamples];
    }

protected:
    std::vector<atlas::math::Point> mSamples;
    atlas::math::Point disk_samples[4];
    std::vector<int> mShuffledIndeces;

    int mNumSamples;
    int mNumSets;
    unsigned long mCount;
    int mJump;
};

class Shape
{
public:
    Shape() : mColour{ 0, 0, 0 }
    {}
    virtual ~Shape() = default;

    // if t computed is less than the t in sr, it and the color should be
    // updated in sr
    virtual bool hit(atlas::math::Ray<atlas::math::Vector> const& ray,
        ShadeRec& sr) const = 0;

    void setColour(Colour const& col)
    {
        mColour = col;
    }

    Colour getColour() const
    {
        return mColour;
    }

    void setMaterial(std::shared_ptr<Material> const& material)
    {
        mMaterial = material;
    }

    std::shared_ptr<Material> getMaterial() const
    {
        return mMaterial;
    }

    virtual bool shadowHit(atlas::math::Ray<atlas::math::Vector> const& ray, float& tMin) const=0;

protected:
    virtual bool intersectRay(atlas::math::Ray<atlas::math::Vector> const& ray,
        float& tMin) const = 0;

    Colour mColour;
    std::shared_ptr<Material> mMaterial;
};

class BRDF
{
public:
    virtual ~BRDF() = default;

    virtual Colour fn(ShadeRec const& sr,
        atlas::math::Vector const& reflected,
        atlas::math::Vector const& incoming) const = 0;
    virtual Colour rho(ShadeRec const& sr,
        atlas::math::Vector const& reflected) const = 0;
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
    
    virtual atlas::math::Vector getDirection(ShadeRec& sr) = 0;

    Colour L([[maybe_unused]] ShadeRec& sr)
    {
        return mRadiance * mColour;
    }

    void scaleRadiance(float b)
    {
        mRadiance = b;
    }

    void setColour(Colour const& c)
    {
        mColour = c;
    }

    bool getShadow() const
    {
        return mShadow;
    }

    void setShadow(const bool shadow)
    {
        mShadow = shadow;
    }

    virtual bool inShadow(atlas::math::Ray<atlas::math::Vector>& ray, ShadeRec& sr) const=0;

protected:
    Colour mColour;
    float mRadiance;
    bool mShadow = true;
};

// Concrete classes which we can construct and use in our ray tracer

class PointLight : public Light
{
public:
    PointLight() : Light{}
    {}

    PointLight(atlas::math::Point const& location) : Light{}
    {
        setLocation(location);
    }

    void setLocation(atlas::math::Vector const& location)
    {
        _location = location;
    }

    atlas::math::Vector getDirection([[maybe_unused]] ShadeRec& sr)
    {
        return glm::normalize(_location - sr.hitPoint);
    }

    bool inShadow(atlas::math::Ray<atlas::math::Vector>& ray, ShadeRec& sr) const
    {
        float t;
        int num_obj = (int)sr.world->scene.size();
        float d = sqrt((_location.x - ray.o.x) * (_location.x - ray.o.x) + (_location.y - ray.o.y) * (_location.y - ray.o.y) + (_location.z - ray.o.z) * (_location.z - ray.o.z));

        for (int i = 0; i < num_obj; i++)
        {
            if (sr.world->scene[i]->shadowHit(ray, t) && t < d)
                return true;
        }
        return false;

    }

private:
    atlas::math::Point _location;
};

class Sphere : public Shape
{
public:
    Sphere(atlas::math::Point center, float radius) :
        mCentre{ center }, mRadius{ radius }, mRadiusSqr{ radius * radius }
    {}

    bool hit(atlas::math::Ray<atlas::math::Vector> const& ray,
        ShadeRec& sr) const
    {
        atlas::math::Vector tmp = ray.o - mCentre;
        float t{ std::numeric_limits<float>::max() };
        bool intersect{ intersectRay(ray, t) };

        // update ShadeRec info about new closest hit
        if (intersect && t < sr.t)
        {
            sr.normal = (tmp + t * ray.d) / mRadius;
            sr.ray = ray;
            sr.color = mColour;
            sr.t = t;
            sr.material = mMaterial;
        }

        return intersect;
    }
    bool shadowHit(atlas::math::Ray<atlas::math::Vector> const& ray, float& tMin) const
    {
        float t;
        const auto tmp{ ray.o - mCentre };
        const auto a{ glm::dot(ray.d, ray.d) };
        const auto b{ 2.0f * glm::dot(ray.d, tmp) };
        const auto c{ glm::dot(tmp, tmp) - mRadiusSqr };
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
    bool intersectRay(atlas::math::Ray<atlas::math::Vector> const& ray,
        float& tMin) const
    {
        const auto tmp{ ray.o - mCentre };
        const auto a{ glm::dot(ray.d, ray.d) };
        const auto b{ 2.0f * glm::dot(ray.d, tmp) };
        const auto c{ glm::dot(tmp, tmp) - mRadiusSqr };
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

    atlas::math::Point mCentre;
    float mRadius;
    float mRadiusSqr;
};

class Triangle : public Shape
{
public:

    Triangle(atlas::math::Point v0, atlas::math::Point v1, atlas::math::Point v2, atlas::math::Normal normal) :
        v0_{ v0 },
        v1_{ v1 },
        v2_{ v2 },
        normal_{ normal }

    {glm::normalize(normal_); }

    // triangle hit function from "Ray Tracing From The Ground Up" textbook
    bool hit(atlas::math::Ray<atlas::math::Vector> const& ray, ShadeRec& sr) const
    {
        float t{ std::numeric_limits<float>::max() };
        bool intersect{ intersectRay(ray, t) };

        // update ShadeRec info about new closest hit
        if (intersect && t < sr.t)
        {
            sr.normal = normal_;
            sr.ray = ray;
            sr.color = mColour;
            sr.t = t;
            sr.material = mMaterial;
        }

        return intersect;
    }

    bool shadowHit(atlas::math::Ray<atlas::math::Vector> const& ray, float& tMin) const
    {
        double a = v0_.x - v1_.x, b = v0_.x - v2_.x, c = ray.d.x, d = v0_.x - ray.o.x;
        double e = v0_.y - v1_.y, f = v0_.y - v2_.y, g = ray.d.y, h = v0_.y - ray.o.y;
        double i = v0_.z - v1_.z, j = v0_.z - v2_.z, k = ray.d.z, l = v0_.z - ray.o.z;

        double m = f * k - g * j, n = h * k - g * l, p = f * l - h * j;
        double q = g * i - e * k, s = e * j - f * i;

        double inv_denom = 1.0 / (a * m + b * q + c * s);

        double e1 = d * m - b * n - c * p;
        double beta = e1 * inv_denom;
        const float kEpsilon{ 0.01f };
        if (beta < kEpsilon)
            return (false);

        double r = e * l - h * i;
        double e2 = a * n + d * q + c * r;
        double gamma = e2 * inv_denom;

        if (gamma < kEpsilon)
            return (false);

        if (beta + gamma > 1.0)
            return (false);

        double e3 = a * p - b * r + d * s;
        double t = e3 * inv_denom;

        if (t < kEpsilon)
            return (false);

        tMin = (float)t;
        return (true);
    }

private:
    bool intersectRay(atlas::math::Ray<atlas::math::Vector> const& ray,
        float& tMin) const
    {
        double a = v0_.x - v1_.x, b = v0_.x - v2_.x, c = ray.d.x, d = v0_.x - ray.o.x;
        double e = v0_.y - v1_.y, f = v0_.y - v2_.y, g = ray.d.y, h = v0_.y - ray.o.y;
        double i = v0_.z - v1_.z, j = v0_.z - v2_.z, k = ray.d.z, l = v0_.z - ray.o.z;

        double m = f * k - g * j, n = h * k - g * l, p = f * l - h * j;
        double q = g * i - e * k, s = e * j - f * i;

        double inv_denom = 1.0 / (a * m + b * q + c * s);

        double e1 = d * m - b * n - c * p;
        double beta = e1 * inv_denom;

        if (beta < 0.0f)
            return (false);

        double r = e * l - h * i;
        double e2 = a * n + d * q + c * r;
        double gamma = e2 * inv_denom;

        if (gamma < 0.0f)
            return (false);

        if (beta + gamma > 1.0f)
            return (false);

        double e3 = a * p - b * r + d * s;
        double t = e3 * inv_denom;

        if (t < 0.0f)
            return (false);

        tMin = (float)t;
        return (true);
    }

    atlas::math::Point v0_, v1_, v2_;
    atlas::math::Normal normal_;
};

class Plane : public Shape
{
public:

    Plane(atlas::math::Point point, atlas::math::Normal normal) :
        point_{ point },
        normal_{ normal }
    {glm::normalize(normal_); }

    bool hit(atlas::math::Ray<atlas::math::Vector> const& ray, ShadeRec& sr) const
    {
        float t{ std::numeric_limits<float>::max() };
        bool intersect{ intersectRay(ray, t) };

        // update ShadeRec info about new closest hit
        if (intersect && t < sr.t)
        {
            sr.normal = normal_;
            sr.ray = ray;
            sr.color = mColour;
            sr.t = t;
            sr.material = mMaterial;
        }

        return intersect;
    }

    bool shadowHit(atlas::math::Ray<atlas::math::Vector> const& ray, float& tMin) const
    {
        float t = glm::dot(point_ - ray.o, normal_) / glm::dot(ray.d, normal_);
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
    bool intersectRay(atlas::math::Ray<atlas::math::Vector> const& ray,
        float& tMin) const
    {
        auto t = glm::dot((point_ - ray.o), normal_) / glm::dot(ray.d, normal_);

        const float kEpsilon{ 0.01f };
        if (atlas::core::geq(t, kEpsilon))
        {
            tMin = t;

            return true;
        }

        return false;
    }

    atlas::math::Point point_;
    atlas::math::Normal normal_;
};

class Regular : public Sampler
{
public:
    Regular(int numSamples, int numSets) : Sampler{ numSamples, numSets }
    {
        generateSamples();
    }

    void generateSamples()
    {
        int n = static_cast<int>(glm::sqrt(static_cast<float>(mNumSamples)));

        for (int j = 0; j < mNumSets; ++j)
        {
            for (int p = 0; p < n; ++p)
            {
                for (int q = 0; q < n; ++q)
                {
                    mSamples.push_back(
                        atlas::math::Point{ (q + 0.5f) / n, (p + 0.5f) / n, 0.0f });
                }
            }
        }
    }
};

class Random : public Sampler
{
public:
    Random(int numSamples, int numSets) : Sampler{ numSamples, numSets }
    {
        generateSamples();
    }

    void generateSamples()
    {
        atlas::math::Random<float> engine;
        for (int p = 0; p < mNumSets; ++p)
        {
            for (int q = 0; q < mNumSamples; ++q)
            {
                mSamples.push_back(atlas::math::Point{
                    engine.getRandomOne(), engine.getRandomOne(), 0.0f });
            }
        }
    }
};

class Lambertian : public BRDF
{
public:
    Lambertian() : mDiffuseColour{}, mDiffuseReflection{}
    {}
    Lambertian(Colour diffuseColor, float diffuseReflection) :
        mDiffuseColour{ diffuseColor }, mDiffuseReflection{ diffuseReflection }
    {}

    Colour
        fn([[maybe_unused]] ShadeRec const& sr,
            [[maybe_unused]] atlas::math::Vector const& reflected,
            [[maybe_unused]] atlas::math::Vector const& incoming) const
    {
        return mDiffuseColour * mDiffuseReflection * glm::one_over_pi<float>();
    }

    Colour
        rho([[maybe_unused]] ShadeRec const& sr,
            [[maybe_unused]] atlas::math::Vector const& reflected) const
    {
        return mDiffuseColour * mDiffuseReflection;
    }

    void setDiffuseReflection(float kd)
    {
        mDiffuseReflection = kd;
    }

    void setDiffuseColour(Colour const& colour)
    {
        mDiffuseColour = colour;
    }

private:
    Colour mDiffuseColour;
    float mDiffuseReflection;
};

class GlossySpecular : public BRDF
{
public:
    GlossySpecular() :ks{}, E{}, cs{}
    {}
    GlossySpecular(float ks_, float exp_, Colour colour):
        ks{ ks_ }, E{ exp_ }, cs{ colour }
    {}

    void set_ks(float ks_)
    {
        ks = ks_;
    }

    void set_exp(float exp_)
    {
        E = exp_;
    }

    void set_cs(Colour const& colour)
    {
        cs = colour;
    }

    Colour
        fn([[maybe_unused]] ShadeRec const& sr,
            [[maybe_unused]] atlas::math::Vector const& wo,
            [[maybe_unused]] atlas::math::Vector const& wi) const
    {
        Colour L{};
        float ndotwi = glm::dot(sr.normal, wi);
        atlas::math::Vector r(-wi + 2.0f * ndotwi * sr.normal);
        
        float rdotwo = glm::dot(r, wo);
        if (rdotwo > 0.0f) {
            L = cs * ks * pow(rdotwo, E);
        }

        return L;
    }

    Colour
        rho([[maybe_unused]] ShadeRec const& sr,
            [[maybe_unused]] atlas::math::Vector const& reflected) const
    {
        return Colour{ 0,0,0 };
    }

private:
    float ks;
    float E;
    Colour cs;

};

class Matte : public Material
{
public:
    Matte() :
        Material{},
        mDiffuseBRDF{ std::make_shared<Lambertian>() },
        mAmbientBRDF{ std::make_shared<Lambertian>() }
    {}

    Matte(float kd, float ka, Colour color) : Matte{}
    {
        setDiffuseReflection(kd);
        setAmbientReflection(ka);
        setDiffuseColour(color);
    }

    void setDiffuseReflection(float k)
    {
        mDiffuseBRDF->setDiffuseReflection(k);
    }

    void setAmbientReflection(float k)
    {
        mAmbientBRDF->setDiffuseReflection(k);
    }

    void setDiffuseColour(Colour colour)
    {
        mDiffuseBRDF->setDiffuseColour(colour);
        mAmbientBRDF->setDiffuseColour(colour);
    }

    Colour shade(ShadeRec& sr)
    {
        using atlas::math::Ray;
        using atlas::math::Vector;

        Vector wo = -sr.ray.o;
        Colour L = mAmbientBRDF->rho(sr, wo) * sr.world->ambient->L(sr);
        size_t numLights = sr.world->lights.size();

        for (size_t i{ 0 }; i < numLights; ++i)
        {
            Vector wi = sr.world->lights[i]->getDirection(sr);
            float nDotWi = glm::dot(sr.normal, wi);

            if (nDotWi > 0.0f)
            {
                bool inShadow = false;
                if (sr.world->lights[i]->getShadow())
                {
                    atlas::math::Ray<atlas::math::Vector> shadowRay(sr.hitPoint, wi);
                    inShadow = sr.world->lights[i]->inShadow(shadowRay, sr);
                }


                
                if (!inShadow) 
                {
                    L += mDiffuseBRDF->fn(sr, wo, wi) * sr.world->lights[i]->L(sr) *
                        nDotWi;
                }
                
            }
        }

        return L;
    }

    

private:
    std::shared_ptr<Lambertian> mDiffuseBRDF;
    std::shared_ptr<Lambertian> mAmbientBRDF;
};

class Phong : public Material
{
public:
    Phong() :
        Material{},
        ambient_brdf{ std::make_shared<Lambertian>() },
        diffuse_brff{ std::make_shared<Lambertian>() },
        specular_brdf{ std::make_shared<GlossySpecular>() }
    {}
    Phong(float ka_,float kd_,  float ks_, float exp_, Colour color) : Phong{}
    {
        setAmbientReflection(ka_); //ka
        setDiffuseReflection(kd_); //kd
        setKsReflection(ks_);
        setExpReflection(exp_);
        setDiffuseColour(color);
    }

    void setAmbientReflection(const float k)
    {
        ambient_brdf->setDiffuseReflection(k);
    }

    void setDiffuseReflection(const float k)
    {
        diffuse_brff->setDiffuseReflection(k);
    }

    void setKsReflection(const float k)
    {
        specular_brdf->set_ks(k);
    }

    void setExpReflection(const float k)
    {
        specular_brdf->set_exp(k);
    }

    void setDiffuseColour(Colour const& colour)
    {
        ambient_brdf->setDiffuseColour(colour);
        diffuse_brff->setDiffuseColour(colour);
        
        specular_brdf->set_cs(colour);
    }

    Colour shade(ShadeRec& sr)
    {
        using atlas::math::Ray;
        using atlas::math::Vector;

        Vector wo = -sr.ray.o;
        Colour L = ambient_brdf->rho(sr, wo) * sr.world->ambient->L(sr);
        size_t numLights = sr.world->lights.size();

        for (size_t i{ 0 }; i < numLights; ++i)
        {
            Vector wi = sr.world->lights[i]->getDirection(sr);
            float nDotWi = glm::dot(sr.normal, wi);

            if (nDotWi > 0.0f)
            {
                bool inShadow = false;
                if (sr.world->lights[i]->getShadow())
                {
                    atlas::math::Ray<atlas::math::Vector> shadowRay(sr.hitPoint, wi);
                    inShadow = sr.world->lights[i]->inShadow(shadowRay, sr);
                }



                if (!inShadow)
                {
                    L += (diffuse_brff->fn(sr, wo, wi) + specular_brdf->fn(sr, wo, wi)) * sr.world->lights[i]->L(sr) * nDotWi;
                }

            }
        }

        return L;
    }

private:
    std::shared_ptr<Lambertian> ambient_brdf;
    std::shared_ptr<Lambertian> diffuse_brff;
    std::shared_ptr<GlossySpecular> specular_brdf;

};

class Directional : public Light
{
public:
    Directional() : Light{}
    {}

    Directional(atlas::math::Vector const& d) : Light{}
    {
        setDirection(d);
    }

    void setDirection(atlas::math::Vector const& d)
    {
        mDirection = glm::normalize(d);
    }

    atlas::math::Vector getDirection([[maybe_unused]] ShadeRec& sr)
    {
        return mDirection;
    }

    bool inShadow([[maybe_unused]] atlas::math::Ray<atlas::math::Vector>& ray, [[maybe_unused]] ShadeRec& sr) const
    {
        return false;
    }

private:
    atlas::math::Vector mDirection;
};

class Ambient : public Light
{
public:
    Ambient() : Light{}
    {}

    atlas::math::Vector getDirection([[maybe_unused]] ShadeRec& sr)
    {
        return atlas::math::Vector{ 0.0f };
    }

    bool inShadow([[maybe_unused]] atlas::math::Ray<atlas::math::Vector>& ray, [[maybe_unused]] ShadeRec& sr) const
    {
        return false;
    }

private:
    atlas::math::Vector mDirection;
};

class Camera
{
public:
    Camera() :
        mEye{ 0.0f, 0.0f, 500.0f },
        mLookAt{ 0.0f },
        mUp{ 0.0f, 1.0f, 0.0f },
        mU{ 1.0f, 0.0f, 0.0f },
        mV{ 0.0f, 1.0f, 0.0f },
        mW{ 0.0f, 0.0f, 1.0f }
    {}

    virtual ~Camera() = default;

    virtual void renderScene(std::shared_ptr<World>& world) const = 0;

    void setEye(atlas::math::Point const& eye)
    {
        mEye = eye;
    }

    void setLookAt(atlas::math::Point const& lookAt)
    {
        mLookAt = lookAt;
    }

    void setUpVector(atlas::math::Vector const& up)
    {
        mUp = up;
    }

    void computeUVW()
    {
        mW = glm::normalize(mEye - mLookAt);
        mU = glm::normalize(glm::cross(mUp, mW));
        mV = glm::cross(mW, mU);

        if (areEqual(mEye.x, mLookAt.x) && areEqual(mEye.z, mLookAt.z) &&
            mEye.y > mLookAt.y)
        {
            mU = { 0.0f, 0.0f, 1.0f };
            mV = { 1.0f, 0.0f, 0.0f };
            mW = { 0.0f, 1.0f, 0.0f };
        }

        if (areEqual(mEye.x, mLookAt.x) && areEqual(mEye.z, mLookAt.z) &&
            mEye.y < mLookAt.y)
        {
            mU = { 1.0f, 0.0f, 0.0f };
            mV = { 0.0f, 0.0f, 1.0f };
            mW = { 0.0f, -1.0f, 0.0f };
        }
    }

protected:
    atlas::math::Point mEye;
    atlas::math::Point mLookAt;
    atlas::math::Point mUp;
    atlas::math::Vector mU, mV, mW;
};

class Pinhole : public Camera
{
public:
    Pinhole::Pinhole() : Camera{}, mDistance{ 500.0f }, mZoom{ 1.0f }
    {}

    void Pinhole::setDistance(float distance)
    {
        mDistance = distance;
    }

    void Pinhole::setZoom(float zoom)
    {
        mZoom = zoom;
    }

    atlas::math::Vector Pinhole::rayDirection(atlas::math::Point const& p) const
    {
        const auto dir = p.x * mU + p.y * mV - mDistance * mW;
        return glm::normalize(dir);
    }

    void Pinhole::renderScene(std::shared_ptr<World>& world) const
    {
        using atlas::math::Point;
        using atlas::math::Ray;
        using atlas::math::Vector;

        Point samplePoint{}, pixelPoint{};
        Ray<atlas::math::Vector> ray{};

        ray.o = mEye;
        float avg{ 1.0f / world->sampler->getNumSamples() };

        for (int r{ 0 }; r < world->height; ++r)
        {
            for (int c{ 0 }; c < world->width; ++c)
            {
                Colour pixelAverage{ 0, 0, 0 };

                for (int j = 0; j < world->sampler->getNumSamples(); ++j)
                {
                    ShadeRec trace_data{};
                    trace_data.world = world;
                    trace_data.t = std::numeric_limits<float>::max();
                    samplePoint = world->sampler->sampleUnitSquare();
                    pixelPoint.x = c - 0.5f * world->width + samplePoint.x;
                    pixelPoint.y = r - 0.5f * world->height + samplePoint.y;
                    ray.d = rayDirection(pixelPoint);

                    bool hit{};
                    for (auto obj : world->scene)
                    {
                        hit |= obj->hit(ray, trace_data);
                    }

                    if (hit)
                    {
                        trace_data.hitPoint = ray.o + trace_data.t * ray.d;
                        pixelAverage += trace_data.material->shade(trace_data);
                    }
                }
                
                pixelAverage.r *= avg;
                pixelAverage.g *= avg;
                pixelAverage.b *= avg;
                // out-of-gamut handling Max-to-one
                if (pixelAverage.r > 1.0f || pixelAverage.g > 1.0f || pixelAverage.b > 1.0f)
                {
                    float MAX_COLOUR = 0.0f;
                    if (pixelAverage.r > MAX_COLOUR)MAX_COLOUR = pixelAverage.r;
                    if (pixelAverage.g > MAX_COLOUR)MAX_COLOUR = pixelAverage.g;
                    if (pixelAverage.b > MAX_COLOUR)MAX_COLOUR = pixelAverage.b;
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
    float mDistance;
    float mZoom;
};

class Fisheye : public Camera
{
public:
    Fisheye() : Camera(),psi_max{60.0f}
    {}

    atlas::math::Vector rayDirection(atlas::math::Point const& p, float& r_squared) const
    {
        atlas::math::Point pn{};
        pn.x = 2.0f / 600 * p.x;
        pn.y = 2.0f / 600 * p.y;
        r_squared = pn.x * pn.x + pn.y * pn.y;

        if (r_squared <= 1.0) {
            float r = sqrt(r_squared);
            float psi = (float)(r * psi_max * 0.017453);
            float sin_psi = sin(psi);
            float cos_psi = cos(psi);
            float sin_alpha = pn.y / r;
            float cos_alpha = pn.x / r;
            const auto dir = sin_psi * cos_alpha * mU + sin_psi * sin_alpha * mV - cos_psi * mW;

            return dir;
        }
        else
            return atlas::math::Vector{0,0,0};
    }

    void renderScene(std::shared_ptr<World>& world) const
    {
        using atlas::math::Point;
        using atlas::math::Ray;
        using atlas::math::Vector;

        Point samplePoint{}, pixelPoint{};
        Ray<atlas::math::Vector> ray{};
        float r_squared;

        ray.o = mEye;
        float avg{ 1.0f / world->sampler->getNumSamples() };

        for (int r{ 0 }; r < world->height; ++r)
        {
            for (int c{ 0 }; c < world->width; ++c)
            {
                Colour pixelAverage{ 0, 0, 0 };

                for (int j = 0; j < world->sampler->getNumSamples(); ++j)
                {
                    ShadeRec trace_data{};
                    trace_data.world = world;
                    trace_data.t = std::numeric_limits<float>::max();
                    samplePoint = world->sampler->sampleUnitSquare();
                    pixelPoint.x = c - 0.5f * world->width + samplePoint.x;
                    pixelPoint.y = r - 0.5f * world->height + samplePoint.y;
                    ray.d = rayDirection(pixelPoint, r_squared);

                    if (r_squared <= 1.0f)
                    {
                        bool hit{};
                        for (auto obj : world->scene)
                        {
                            hit |= obj->hit(ray, trace_data);
                        }

                        if (hit)
                        {
                            trace_data.hitPoint = ray.o + trace_data.t * ray.d;
                            pixelAverage += trace_data.material->shade(trace_data);
                        }
                    }
                }

                pixelAverage.r *= avg;
                pixelAverage.g *= avg;
                pixelAverage.b *= avg;
                // out-of-gamut handling Max-to-one
                if (pixelAverage.r > 1.0f || pixelAverage.g > 1.0f || pixelAverage.b > 1.0f)
                {
                    float MAX_COLOUR = 0.0f;
                    if (pixelAverage.r > MAX_COLOUR)MAX_COLOUR = pixelAverage.r;
                    if (pixelAverage.g > MAX_COLOUR)MAX_COLOUR = pixelAverage.g;
                    if (pixelAverage.b > MAX_COLOUR)MAX_COLOUR = pixelAverage.b;
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
    float psi_max;

};

class ThinLens : public Camera
{
public:
    ThinLens() : Camera{},lens_radius{ 0.25f }, d{ 30.0f }, f{ 25.0f }, zoom{ 1.0f }, sampler_ptr{std::make_shared<Regular>(4, 1) }
    {sampler_ptr->map_samples_to_unit_disk(); }


    atlas::math::Vector rayDirection(atlas::math::Point const& pixel_point, atlas::math::Point const& lens_point) const
    {
        atlas::math::Point p{};
        p.x = pixel_point.x * f / d;
        p.y = pixel_point.y * f / d;

        atlas::math::Vector dir = (p.x - lens_point.x) * mU + (p.y - lens_point.y) * mV - f * mW;
        return glm::normalize(dir);
    }

    void renderScene(std::shared_ptr<World>& world) const
    {
        using atlas::math::Point;
        using atlas::math::Ray;
        using atlas::math::Vector;

        Point sp{};
        Point pp{};
        Point dp{};
        Point lp{};
        Ray<atlas::math::Vector> ray{};

        //ray.o = mEye;
        float avg{ 1.0f / world->sampler->getNumSamples() };

        for (int r{ 0 }; r < world->height; ++r)
        {
            for (int c{ 0 }; c < world->width; ++c)
            {
                Colour pixelAverage{ 0, 0, 0 };

                for (int j = 0; j < world->sampler->getNumSamples(); ++j)
                {
                    ShadeRec trace_data{};
                    trace_data.world = world;
                    trace_data.t = std::numeric_limits<float>::max();
                    sp = world->sampler->sampleUnitSquare();
                    pp.x = c - 0.5f * world->width + sp.x;
                    pp.y = r - 0.5f * world->height + sp.y;

                    dp = sampler_ptr->sampleUnitDisk();
                    lp = dp * lens_radius;

                    ray.o = mEye + lp.x * mU + lp.y * mV;
                    ray.d = rayDirection(pp, lp);

                    bool hit{};
                    for (auto obj : world->scene)
                    {
                        hit |= obj->hit(ray, trace_data);
                    }

                    if (hit)
                    {
                        trace_data.hitPoint = ray.o + trace_data.t * ray.d;
                        pixelAverage += trace_data.material->shade(trace_data);
                    }
                }

                pixelAverage.r *= avg;
                pixelAverage.g *= avg;
                pixelAverage.b *= avg;
                // out-of-gamut handling Max-to-one
                if (pixelAverage.r > 1.0f || pixelAverage.g > 1.0f || pixelAverage.b > 1.0f)
                {
                    float MAX_COLOUR = 0.0f;
                    if (pixelAverage.r > MAX_COLOUR)MAX_COLOUR = pixelAverage.r;
                    if (pixelAverage.g > MAX_COLOUR)MAX_COLOUR = pixelAverage.g;
                    if (pixelAverage.b > MAX_COLOUR)MAX_COLOUR = pixelAverage.b;
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
    float lens_radius;
    float d;
    float f;
    float zoom;
    std::shared_ptr<Sampler> sampler_ptr;
};