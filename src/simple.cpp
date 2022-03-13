#include "nori/color.h"
#include "nori/common.h"
#include "nori/mesh.h"
#include "nori/sampler.h"
#include "nori/vector.h"
#include <nori/integrator.h>
#include <nori/scene.h>
NORI_NAMESPACE_BEGIN

class SimpleIntegrator : public Integrator {
    Point3f lightPos;
    Color3f lightEnergy;
public:
    SimpleIntegrator(const PropertyList& props)
    {
        lightPos = props.getPoint("position");
        lightEnergy = props.getColor("energy");
    }
    Color3f Li(const Scene* scene, Sampler* sampler, const Ray3f& ray) const
    {
        
        Intersection its;
        if (!scene->rayIntersect(ray, its)) {
            return Color3f(0.0);
        } else {
            Point3f shadingPoint = its.p;
            Normal3f normal = its.shFrame.n;
            Vector3f outDir = lightPos - shadingPoint;
            float squareOfDir = outDir.dot(outDir);
            Color3f result;
            Intersection shadowRayIts;
            if (scene->shadowrayIntersect(Ray3f(shadingPoint, outDir),shadowRayIts)) {
                result = Color3f(0, 0, 0);
            } else {
                result = .25 * INV_PI * INV_PI * lightEnergy * fmaxf(0.0, normal.dot(outDir.normalized())) / squareOfDir;
            }
            return result;
        }
    }
    /// Return a human-readable description for debugging purposes
    std::string toString() const
    {
        return "SimpleIntegrator[]";
    }
};
NORI_REGISTER_CLASS(SimpleIntegrator, "simple");
NORI_NAMESPACE_END