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
            Point3f shadingPoint = its.p;
            Normal3f shadingPointNormal = its.shFrame.n;
            const BSDF* bsdf = its.mesh->getBSDF();
            Vector3f wiLocal = its.shFrame.toLocal(-ray.d);
            BSDFQueryRecord bsdfQueryRecord { wiLocal };
            Color3f result { 0., 0., 0. };
            if (its.mesh->isEmitter()) {
                result = its.mesh->getEmitter()->getRadiance();
            } else if (bsdf->isDiffuse()) {
                float emitterPdf;
                auto mesh = scene->sampleEmitter(sampler, emitterPdf);
                if (mesh != nullptr) {
                    Point3f lightSamplePos;
                    Vector3f lightSamplePosNormal;
                    float pdf;
                    mesh->sample(sampler, lightSamplePos, lightSamplePosNormal, pdf);
                    Vector3f outDir = lightSamplePos - shadingPoint;
                    float distanceSquared = outDir.dot(outDir);
                    //算完距离之后单位化
                    outDir.normalize();
                    bsdfQueryRecord.wo = its.shFrame.toLocal(outDir);
                    bsdfQueryRecord.measure = ESolidAngle;
                    Color3f bsdfValue = bsdf->eval(bsdfQueryRecord);
                    auto secondaryRay = Ray3f(shadingPoint, outDir);
                    Intersection shadowRayIts;
                    scene->shadowrayIntersect(secondaryRay, shadowRayIts);
                    if (shadowRayIts.t * shadowRayIts.t >= distanceSquared - Epsilon) {
                        result = bsdfValue * mesh->getEmitter()->sample(-outDir, lightSamplePosNormal, distanceSquared) * fmaxf(0.0, shadingPointNormal.dot(outDir)) / (pdf * emitterPdf);
                    }
                }
            } else {
                const float russianRouletteProbability = 0.95;
                if (sampler->next1D() < russianRouletteProbability) {
                    bsdf->sample(bsdfQueryRecord, sampler->next2D());
                    Vector3f woWorld = its.shFrame.toWorld(bsdfQueryRecord.wo);
                    result = Li(scene, sampler, Ray3f { shadingPoint, woWorld }) / russianRouletteProbability;
                }
            }
            return result;
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