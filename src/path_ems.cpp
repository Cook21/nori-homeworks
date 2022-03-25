#include "nori/color.h"
#include "nori/common.h"
#include "nori/dpdf.h"
#include "nori/mesh.h"
#include "nori/sampler.h"
#include "nori/vector.h"
#include <functional>
#include <nori/bsdf.h>
#include <nori/emitter.h>
#include <nori/integrator.h>
#include <nori/scene.h>
NORI_NAMESPACE_BEGIN
///Material Path Tracer
class PathEmsIntegrator : public Integrator {

public:
    PathEmsIntegrator(const PropertyList& props)
    {
    }
    Color3f Li(const Scene* scene, Sampler* sampler, const Ray3f& ray) const
    {
        return LiForNEE(scene, sampler, ray, false);
    }
    Color3f LiForNEE(const Scene* scene, Sampler* sampler, const Ray3f& ray, bool isDiffuseBounce) const
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
            if (its.mesh->isEmitter() && !isDiffuseBounce) { 
                result += its.mesh->getEmitter()->getRadiance();
            }
            
            //计算间接光和反射光
            float russianRouletteProbability;
            if (bsdf->isDiffuse()) {
                russianRouletteProbability = 0.75;
            } else {
                russianRouletteProbability = 0.95;
            }
            if (sampler->next1D() < russianRouletteProbability) {
                Color3f bsdfResult = bsdf->sample(bsdfQueryRecord, sampler->next2D());
                Vector3f woWorld = its.shFrame.toWorld(bsdfQueryRecord.wo);
                result += bsdfResult * LiForNEE(scene, sampler, Ray3f { shadingPoint, woWorld }, bsdf->isDiffuse()) / russianRouletteProbability;
            }
            //计算直接光
            if (bsdf->isDiffuse()) {
                EmitterQueryRecord eRec = EmitterQueryRecord(ESolidAngle);
                float pdf;
                result += scene->sampleEmitter(sampler, bsdf, bsdfQueryRecord, its, eRec, pdf);
            }
            return result;
        }
    }
    /// Return a human-readable description for debugging purposes
    std::string toString() const
    {
        return "PathEmsIntegrator[]";
    }
};
NORI_REGISTER_CLASS(PathEmsIntegrator, "path_ems");
NORI_NAMESPACE_END