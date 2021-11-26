#include "nori/color.h"
#include "nori/common.h"
#include "nori/mesh.h"
#include "nori/sampler.h"
#include "nori/vector.h"
#include <nori/scene.h>
#include <nori/integrator.h>
NORI_NAMESPACE_BEGIN

class NormalIntegrator : public Integrator{
    public:
    NormalIntegrator(const PropertyList &props){
        
    }
    Color3f Li(const Scene* scene, Sampler* sampler,const Ray3f& ray) const {
        Intersection its;
        if(!scene->rayIntersect(ray,its)){
            return Color3f(0.0);
        }else{
            Normal3f normal=its.shFrame.n.cwiseAbs();
            return Color3f(normal.x(),normal.y(),normal.z());
        }
        
    }
    /// Return a human-readable description for debugging purposes
    std::string toString() const {
        return "NormalIntegrator[]";
    }
};
NORI_REGISTER_CLASS(NormalIntegrator, "normals");
NORI_NAMESPACE_END