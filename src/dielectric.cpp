/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob

    Nori is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Nori is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <nori/bsdf.h>
#include <nori/common.h>
#include <nori/frame.h>

NORI_NAMESPACE_BEGIN

/// Ideal dielectric BSDF
class Dielectric : public BSDF {
public:
    Dielectric(const PropertyList& propList)
    {
        //IOR = Index Of Refraction
        /* Interior IOR (default: BK7 borosilicate optical glass) */
        m_intIOR = propList.getFloat("intIOR", 1.5046f);

        /* Exterior IOR (default: air) */
        m_extIOR = propList.getFloat("extIOR", 1.000277f);
    }

    Color3f eval(const BSDFQueryRecord&) const
    {
        /* Discrete BRDFs always evaluate to zero in Nori */
        return Color3f(0.0f);
    }

    float pdf(const BSDFQueryRecord&) const
    {
        /* Discrete BRDFs always evaluate to zero in Nori */
        return 0.0f;
    }

    Color3f sample(BSDFQueryRecord& bRec, const Point2f& sample) const
    {
        //throw NoriException("Unimplemented!");
        Vector3f normalLocal { 0., 0., 1. };
        float cosThetaI = Frame::cosTheta(bRec.wi);
        bool entering = (cosThetaI > 0.f);
        float etaI, etaT;
        if (entering) {
            etaI = m_extIOR;
            etaT = m_intIOR;
        } else {
            etaI = m_intIOR;
            etaT = m_extIOR;
            normalLocal = -normalLocal;
            cosThetaI = -cosThetaI;
        }
        bRec.eta = etaI / etaT;
        float sin2ThetaI = Frame::sinTheta2(bRec.wi);
        float sin2ThetaT = bRec.eta * bRec.eta * sin2ThetaI;
        if (sin2ThetaT >= 1.) { //全反射
            bRec.wo = Vector3f { -bRec.wi.x(), -bRec.wi.y(), bRec.wi.z() };
        } else {
            //反射的比例，精确
            float fresnelTerm = nori::fresnel(cosThetaI, etaI, etaT);
            //反射的比例，近似
            //float fre =  SchlicksFresnel(cosThetaI, etaI, etaT);
            if (sample.x() < fresnelTerm) {
                bRec.wo = Vector3f { -bRec.wi.x(), -bRec.wi.y(), bRec.wi.z() };
            } else {
                float cosThetaT = std::sqrt(1. - sin2ThetaT);
                //球坐标方法
                
                float thetaT;
                if (entering) {
                    thetaT = M_PI - acos(cosThetaT);
                } else {
                    thetaT = acos(cosThetaT);
                }
                float phiT = std::atan2(bRec.wi.y(), bRec.wi.x()) + M_PI;
                bRec.wo = sphericalDirection(thetaT, phiT);
                
                //向量方法
                //bRec.wo = bRec.eta * -bRec.wi + (bRec.eta * cosThetaI - cosThetaT) * normalLocal;
            }
        }
        return Color3f(1.);
    }

    std::string toString() const
    {
        return tfm::format(
            "Dielectric[\n"
            "  intIOR = %f,\n"
            "  extIOR = %f\n"
            "]",
            m_intIOR, m_extIOR);
    }

private:
    float m_intIOR, m_extIOR;
    //Schlick’s approximation
    float SchlicksFresnel(float cosThetaI, float etaI, float etaT) const
    {
        float temp = (etaI - etaT) / (etaI + etaT);
        float r0 = temp * temp;
        return r0 + (1. - r0) * pow((1 - cosThetaI), 5.);
    }
};

NORI_REGISTER_CLASS(Dielectric, "dielectric");
NORI_NAMESPACE_END
