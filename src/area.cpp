#include "nori/color.h"
#include "nori/emitter.h"
#include "nori/common.h"
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
    Color3f sample(Vector3f wiWorld, Vector3f normalWorld, float distanceSquared) const override {
        Color3f result = radiance * std::max(wiWorld.dot(normalWorld),0.0f) / distanceSquared;
        return result;
    }
    Color3f getRadiance() const {
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