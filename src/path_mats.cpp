#include "nori/color.h"
#include "nori/common.h"
#include "nori/mesh.h"
#include "nori/sampler.h"
#include "nori/vector.h"
#include <nori/bsdf.h>
#include <nori/emitter.h>
#include <nori/integrator.h>
#include <nori/scene.h>
NORI_NAMESPACE_BEGIN
///Material Path Tracer
class PathMatsIntegrator : public Integrator {

public:
    PathMatsIntegrator(const PropertyList& props)
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
            const BSDF* bsdf = its.mesh->getBSDF();
            Vector3f wiLocal = its.shFrame.toLocal(-ray.d);
            BSDFQueryRecord bsdfQueryRecord { wiLocal };
            Color3f result { 0., 0., 0. };
            if (its.mesh->isEmitter()) {
                result += its.mesh->getEmitter()->getRadiance();
            }
            float russianRouletteProbability;
            if (bsdf->isDiffuse()) {
                russianRouletteProbability = 0.75;
            } else {
                russianRouletteProbability = 0.95;
            }
            if (sampler->next1D() < russianRouletteProbability) {
                Color3f bsdfResult = bsdf->sample(bsdfQueryRecord, sampler->next2D());
                Vector3f woWorld = its.shFrame.toWorld(bsdfQueryRecord.wo);
                result += bsdfResult * Li(scene, sampler, Ray3f { shadingPoint, woWorld }) / russianRouletteProbability;
            }

            return result;
        }
    }
    /// Return a human-readable description for debugging purposes
    std::string toString() const
    {
        return "PathMatsIntegrator[]";
    }
};
NORI_REGISTER_CLASS(PathMatsIntegrator, "path_mats");
NORI_NAMESPACE_END