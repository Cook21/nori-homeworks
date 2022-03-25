#include "nori/color.h"
#include "nori/common.h"
#include "nori/dpdf.h"
#include "nori/mesh.h"
#include "nori/proplist.h"
#include "nori/sampler.h"
#include "nori/vector.h"
#include <functional>
#include <nori/bsdf.h>
#include <nori/emitter.h>
#include <nori/integrator.h>
#include <nori/scene.h>
NORI_NAMESPACE_BEGIN
///Material Path Tracer
class PathMisIntegrator : public Integrator {

public:
    PathMisIntegrator(const PropertyList& props)
    {
    }
    Color3f Li(const Scene* scene, Sampler* sampler, const Ray3f& ray) const
    {
        return LiForMIS(scene, sampler, ray, NAN);
    }
    Color3f LiForMIS(const Scene* scene, Sampler* sampler, const Ray3f& ray, float pBSDF) const
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
                float weight = 1.0f;
                if (isnan(pBSDF)) {
                    weight = 1.0f;
                }else{
                    EmitterQueryRecord eRec { ESolidAngle, wiLocal, its.t * its.t };
                    float pLight = scene->getEmitterPdf(eRec, its.mesh);
                    weight = pBSDF / (pLight + pBSDF);
                }
                result += its.mesh->getEmitter()->getRadiance() * weight;
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
                float pBSDF;
                if (bsdf->isDiffuse())
                    pBSDF = bsdf->pdf(bsdfQueryRecord);
                else{
                    pBSDF=NAN;
                }
                result += bsdfResult * LiForMIS(scene, sampler, Ray3f { shadingPoint, woWorld }, pBSDF) / russianRouletteProbability;
            }
            //计算直接光
            if (bsdf->isDiffuse()) {
                EmitterQueryRecord eRec = EmitterQueryRecord(ESolidAngle);
                float pLight;
                Color3f lightSampleResult = scene->sampleEmitter(sampler, bsdf, bsdfQueryRecord, its, eRec, pLight);
                if (isnan(pLight)) {
                    //occluded, do nothing
                } else {
                    float pBSDF = bsdf->pdf(bsdfQueryRecord);
                    float weight;
                    if(isinf(pLight)){
                        weight=1.0f;
                    }else{
                        weight = pLight / (pLight + pBSDF);
                    }
                    result += lightSampleResult * weight;
                }
            }
            return result;
        }
    }
    /// Return a human-readable description for debugging purposes
    std::string toString() const
    {
        return "PathMisIntegrator[]";
    }
};
NORI_REGISTER_CLASS(PathMisIntegrator, "path_mis");
NORI_NAMESPACE_END