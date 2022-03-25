#include "nori/color.h"
#include "nori/common.h"
#include "nori/emitter.h"
#include "nori/frame.h"
#include "nori/mesh.h"
#include "nori/vector.h"
#include <nori/scene.h>
NORI_NAMESPACE_BEGIN
class AreaLight : public Emitter {
    Color3f radiance;

public:
    AreaLight(const PropertyList& props)
    {
        radiance = props.getColor("radiance");
    }
    float pdf(const EmitterQueryRecord& eRec, const float pdfIn) const override{
        if(eRec.targetMeasure == EDiscrete){
            return pdfIn;
        }else if(eRec.targetMeasure == ESolidAngle){
            float areaToSolidAngleFactor = std::max(Frame::cosTheta(eRec.wiLocal), 0.0f) / eRec.distanceSquared;
            return pdfIn / areaToSolidAngleFactor;
        }else{
            return NAN;
        }
    }
    Color3f getRadiance() const override
    {
        return radiance;
    }
    /// Return a human-readable description for debugging purposes
    std::string toString() const override
    {
        return "AreaLight[]";
    }
};
NORI_REGISTER_CLASS(AreaLight, "area");
NORI_NAMESPACE_END