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

#include "nori/common.h"
#include <nori/bsdf.h>
#include <nori/frame.h>
#include <nori/warp.h>

NORI_NAMESPACE_BEGIN

class Microfacet : public BSDF {
public:
    Microfacet(const PropertyList& propList)
    {
        /* RMS surface roughness */
        m_alpha = propList.getFloat("alpha", 0.1f);

        /* Interior IOR (default: BK7 borosilicate optical glass) */
        m_intIOR = propList.getFloat("intIOR", 1.5046f);

        /* Exterior IOR (default: air) */
        m_extIOR = propList.getFloat("extIOR", 1.000277f);

        /* Albedo of the diffuse base material (a.k.a "kd") */
        m_kd = propList.getColor("kd", Color3f(0.5f));

        /* To ensure energy conservation, we must scale the 
           specular component by 1-kd. 

           While that is not a particularly realistic model of what 
           happens in reality, this will greatly simplify the 
           implementation. Please see the course staff if you're 
           interested in implementing a more realistic version 
           of this BRDF. */
        m_ks = 1 - m_kd.maxCoeff();
    }

    /// Evaluate the BRDF for the given pair of directions
    Color3f eval(const BSDFQueryRecord& bRec) const
    {
        if (Frame::cosTheta(bRec.wi) > 0 && Frame::cosTheta(bRec.wo)>0 && bRec.measure == ESolidAngle) {
            auto G1 = [](const float& alpha, Vector3f w1, Vector3f w2) -> float {
                const Vector3f n = Vector3f(0., 0., 1.);
                float c = w1.dot(w2) / w1.dot(n);
                if (c > 0) {
                    auto b = 1 / (alpha * Frame::tanTheta(w1));
                    if (b < 1.6) {
                        return b * (3.535 + 2.181 * b) / (1. + b * (2.276 + 2.577 * b));
                    } else {
                        return  1.;
                    }
                } else {
                    return 0.;
                }
            };
            const Vector3f &wi = bRec.wi, &wo = bRec.wo;
            Vector3f wh = (bRec.wi + bRec.wo).normalized();
            float D = Warp::squareToBeckmannPdf(wh, m_alpha);
            float F = fresnel(wh.dot(wi), m_extIOR, m_intIOR);
            float G = G1(m_alpha, wi, wh) * G1(m_alpha, wo, wh);
            return INV_PI * m_kd + m_ks * D * F * G / (4. * Frame::cosTheta(wi) * Frame::cosTheta(wo) * Frame::cosTheta(wh));
        }else{
            return 0.f;
        }
    }

    /// Evaluate the sampling density of \ref sample() wrt. solid angles
    float pdf(const BSDFQueryRecord& bRec) const
    {
        if (Frame::cosTheta(bRec.wi) > 0 && Frame::cosTheta(bRec.wo)>0 && bRec.measure == ESolidAngle) {
            Vector3f wh = (bRec.wi + bRec.wo).normalized();
            float Jacobian = 1 / (4 * wh.dot(bRec.wo));
            return m_ks * Warp::squareToBeckmannPdf(wh, m_alpha) * Jacobian + INV_PI * (1. - m_ks) * Frame::cosTheta(bRec.wo);
        } else {
            return 0.f;
        }
    }

    /// Sample the BRDF
    Color3f sample(BSDFQueryRecord& bRec, const Point2f& _sample) const
    {
        //throw NoriException("MicrofacetBRDF::sample(): not implemented!");

        // Note: Once you have implemented the part that computes the scattered
        // direction, the last part of this function should simply return the
        // BRDF value divided by the solid angle density and multiplied by the
        // cosine factor from the reflection equation, i.e.
        if(Frame::cosTheta(bRec.wi)<=0){
            return Color3f(.0f);
        }
        else if(_sample.x()<m_ks){ //Specular case
            Point2f sample {_sample.x()/m_ks,_sample.y()}; //sample reuse
            Vector3f normal = Warp::squareToBeckmann(sample,m_alpha);
            bRec.wo = 2.f * normal.dot(bRec.wi) * normal - bRec.wi;
            //检查数值稳定性
            // if((normal - (bRec.wo+bRec.wi).normalized()).norm() > 0.5f)
            //     cout <<"wi:"<< bRec.wi << "normal:"<< normal << "wo" << bRec.wo<< endl;
        }else{ //Diffuse case
            Point2f sample {(_sample.x()-m_ks)/(1.f-m_ks),_sample.y()}; //sample reuse
            bRec.wo = Warp::squareToCosineHemisphere(sample);
        }
        bRec.measure = ESolidAngle;
        float pdf = this->pdf(bRec);
        if(pdf>0)
            return eval(bRec) * Frame::cosTheta(bRec.wo) / pdf;
        else{
            return Color3f(0.f);
        }
    }

    bool isDiffuse() const
    {
        /* While microfacet BRDFs are not perfectly diffuse, they can be
           handled by sampling techniques for diffuse/non-specular materials,
           hence we return true here */
        return true;
    }

    std::string toString() const
    {
        return tfm::format(
            "Microfacet[\n"
            "  alpha = %f,\n"
            "  intIOR = %f,\n"
            "  extIOR = %f,\n"
            "  kd = %s,\n"
            "  ks = %f\n"
            "]",
            m_alpha,
            m_intIOR,
            m_extIOR,
            m_kd.toString(),
            m_ks);
    }

private:
    float m_alpha;
    float m_intIOR, m_extIOR;
    float m_ks;
    Color3f m_kd;
};

NORI_REGISTER_CLASS(Microfacet, "microfacet");
NORI_NAMESPACE_END
