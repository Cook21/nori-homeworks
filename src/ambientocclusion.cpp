#include "nori/color.h"
#include "nori/common.h"
#include "nori/mesh.h"
#include "nori/sampler.h"
#include "nori/vector.h"
#include "nori/warp.h"
#include <nori/integrator.h>
#include <nori/scene.h>
#include <pcg32.h>
NORI_NAMESPACE_BEGIN

class AmbientOcclusionIntegrator : public Integrator {
    Color3f ambientRadiance = { 1.0, 1.0, 1.0 };

public:
    AmbientOcclusionIntegrator(const PropertyList& props)
    {
    }
    Color3f Li(const Scene* scene, Sampler* sampler, const Ray3f& ray) const
    {

        Intersection its;
        if (!scene->rayIntersect(ray, its)) {
            return Color3f(0.0);
        } else {
            Point3f shadingPoint = its.p;
            Normal3f normal = its.shFrame.n;
            // if(abs(normal.norm()-1.f)>.001f)
            //     cout << "normal length incorrected!\n";
            Color3f result { 0, 0, 0 };
            auto seed = sampler->next2D();
            auto sample = Warp::squareToCosineHemisphere(seed);
            auto probabilityDensity = Warp::squareToCosineHemispherePdf(sample);
            Vector3f outDir = its.shFrame.toWorld(sample);
            if (!scene->rayIntersect(Ray3f(shadingPoint, outDir))) {
                //Monte Carlo积分，要除以PDF
                result = INV_PI * ambientRadiance * std::max(normal.dot(outDir.normalized()),0.0f) / probabilityDensity;
            }
            return result;
        }
    }
    /// Return a human-readable description for debugging purposes
    std::string toString() const
    {
        return "AmbientOcclusionIntegrator[]";
    }
};
NORI_REGISTER_CLASS(AmbientOcclusionIntegrator, "ao");
NORI_NAMESPACE_END