#include "nori/color.h"
#include "nori/common.h"
#include "nori/mesh.h"
#include "nori/sampler.h"
#include "nori/vector.h"
#include <nori/emitter.h>
#include <nori/integrator.h>
#include <nori/scene.h>
NORI_NAMESPACE_BEGIN

class WhittedIntegrator : public Integrator {

public:
    WhittedIntegrator(const PropertyList& props)
    {
    }
    Color3f Li(const Scene* scene, Sampler* sampler, const Ray3f& ray) const
    {
        Intersection its;
        if (!scene->rayIntersect(ray, its)) {
            return Color3f(0.0);
        } else {
            Point3f shadingPoint = its.p;
            Normal3f shadingPointNormal = its.shFrame.n;
            Color3f result{0.,0.,0.};
            int count=0;
            for (auto mesh : scene->getMeshes()) {
                if (mesh->isEmitter()) {
                    Point3f lightSamplePos;
                    Vector3f lightSamplePosNormal;
                    float pdf;
                    mesh->sample(sampler, lightSamplePos, lightSamplePosNormal, pdf);
                    Vector3f outDir = lightSamplePos - shadingPoint;
                    float distanceSquared = outDir.dot(outDir);
                    //算完距离之后单位化
                    outDir.normalize();
                    Intersection its;
                    scene->rayIntersect(Ray3f(shadingPoint, outDir), its);
                    if (its.mesh == mesh) {
                        result += mesh->getEmitter()->sample(-outDir, lightSamplePosNormal, distanceSquared) * fmaxf(0.0, shadingPointNormal.dot(outDir)) / pdf ;
                        count++;
                    }
                }
            }
            return result / count;
        }
    }
    /// Return a human-readable description for debugging purposes
    std::string toString() const
    {
        return "WhittedIntegrator[]";
    }
};
NORI_REGISTER_CLASS(WhittedIntegrator, "whitted");
NORI_NAMESPACE_END