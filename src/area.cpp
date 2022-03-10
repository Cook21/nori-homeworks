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
    Color3f sample(Vector3f wiWorld, Vector3f normalWorld, float distanceSquared) override {
        Color3f result = radiance * wiWorld.dot(normalWorld) / distanceSquared;
        return result;
    }
    /// Return a human-readable description for debugging purposes
    std::string toString() const override
    {
        return "AreaLight[]";
    }
};
NORI_REGISTER_CLASS(AreaLight, "area");
NORI_NAMESPACE_END