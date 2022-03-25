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
            Point3f& shadingPoint = its.p;
            const BSDF* bsdf = its.mesh->getBSDF();
            Vector3f wiLocal = its.shFrame.toLocal(-ray.d);
            BSDFQueryRecord bsdfQueryRecord { wiLocal };
            Color3f result { 0., 0., 0. };
            if (its.mesh->isEmitter()) {
                result += its.mesh->getEmitter()->getRadiance();
            }
            if (bsdf->isDiffuse()) {
                EmitterQueryRecord eRec = EmitterQueryRecord(ESolidAngle);
                float pdf;
                result += scene->sampleEmitter(sampler, bsdf, bsdfQueryRecord, its, eRec, pdf);
            } else {
                const float russianRouletteProbability = 0.95;
                if (sampler->next1D() < russianRouletteProbability) {
                    bsdf->sample(bsdfQueryRecord, sampler->next2D());
                    Vector3f woWorld = its.shFrame.toWorld(bsdfQueryRecord.wo);
                    result += Li(scene, sampler, Ray3f { shadingPoint, woWorld }) / russianRouletteProbability;
                }
            }
            return result;
        }
    }
    /// Return a human-readable description for debugging purposes
    std::string
    toString() const
    {
        return "WhittedIntegrator[]";
    }
};
NORI_REGISTER_CLASS(WhittedIntegrator, "whitted");
NORI_NAMESPACE_END