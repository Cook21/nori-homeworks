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

#pragma once

#include "nori/bsdf.h"
#include "nori/color.h"
#include "nori/common.h"
#include "nori/frame.h"
#include "nori/mesh.h"
#include <cstddef>
#include <nori/accel.h>
#include <nori/emitter.h>
#include <nori/sampler.h>

NORI_NAMESPACE_BEGIN

/**
 * \brief Main scene data structure
 *
 * This class holds information on scene objects and is responsible for
 * coordinating rendering jobs. It also provides useful query routines that
 * are mostly used by the \ref Integrator implementations.
 */
class Scene : public NoriObject {
public:
    /// Construct a new scene object
    Scene(const PropertyList&);

    /// Release all memory
    virtual ~Scene();

    /// Return a pointer to the scene's kd-tree
    const Accel* getAccel() const { return m_accel; }

    /// Return a pointer to the scene's integrator
    const Integrator* getIntegrator() const { return m_integrator; }

    /// Return a pointer to the scene's integrator
    Integrator* getIntegrator() { return m_integrator; }

    /// Return a pointer to the scene's camera
    const Camera* getCamera() const { return m_camera; }

    /// Return a pointer to the scene's sample generator (const version)
    const Sampler* getSampler() const { return m_sampler; }

    /// Return a pointer to the scene's sample generator
    Sampler* getSampler() { return m_sampler; }

    /// Return a reference to an array containing all meshes
    const std::vector<Mesh*>& getMeshes() const { return m_meshes; }

    /**
     * \brief Intersect a ray against all triangles stored in the scene
     * and return detailed intersection information
     *
     * \param ray
     *    A 3-dimensional ray data structure with minimum/maximum
     *    extent information
     *
     * \param its
     *    A detailed intersection record, which will be filled by the
     *    intersection query
     *
     * \return \c true if an intersection was found
     */
    bool rayIntersect(const Ray3f& ray, Intersection& its) const
    {
        return m_accel->rayIntersect(ray, its, false);
    }

    /**
     * \brief Intersect a ray against all triangles stored in the scene
     * and \a only determine whether or not there is an intersection.
     *
     * This method much faster than the other ray tracing function,
     * but the performance comes at the cost of not providing any
     * additional information about the detected intersection
     * (not even its position).
     *
     * \param ray
     *    A 3-dimensional ray data structure with minimum/maximum
     *    extent information
     *
     * \return \c true if an intersection was found
     */
    bool shadowrayIntersect(const Ray3f& ray, Intersection& its) const
    {
        return m_accel->rayIntersect(ray, its, true);
    }

    /// \brief Return an axis-aligned box that bounds the scene
    const BoundingBox3f& getBoundingBox() const
    {
        return m_accel->getBoundingBox();
    }

    /**
     * \brief Inherited from \ref NoriObject::activate()
     *
     * Initializes the internal data structures (kd-tree,
     * emitter sampling data structures, etc.)
     */
    void activate();

    /// Add a child object to the scene (meshes, integrators etc.)
    void addChild(NoriObject* obj);

    /// Return a string summary of the scene (for debugging purposes)
    std::string toString() const;

    Color3f sampleEmitter(Sampler* sampler, const BSDF* bsdf, BSDFQueryRecord& bRec, const Intersection& its, EmitterQueryRecord& eRec, float& pdfOut) const
    {
        Color3f result { .0f };
        if (emitterDPDF.size() > 0) {
            float emitterIdPdf;
            size_t emitterId = emitterDPDF.sample(sampler->next1D(), emitterIdPdf);
            Mesh* mesh = m_meshes[emitterMeshId[emitterId]];
            Point3f lightSamplePos;
            float SamplePosPdf;
            //生成法线
            Vector3f normalWorld;
            mesh->sample(sampler, lightSamplePos, normalWorld, SamplePosPdf);
            Vector3f outDir = lightSamplePos - its.p;
            eRec.distanceSquared = outDir.squaredNorm();
            //算完距离之后单位化
            outDir.normalize();
            nori::Frame frame { normalWorld };
            eRec.wiLocal = frame.toLocal(-outDir);
            auto shadowRay = Ray3f(its.p, outDir);
            Intersection shadowRayIts;
            this->shadowrayIntersect(shadowRay, shadowRayIts);
            if (shadowRayIts.t * shadowRayIts.t >= eRec.distanceSquared - Epsilon) {
                bRec.wo = its.shFrame.toLocal(outDir);
                bRec.measure = ESolidAngle;
                Color3f bsdfValue = bsdf->eval(bRec);
                pdfOut = mesh->getEmitter()->pdf(eRec, emitterIdPdf * SamplePosPdf);
                result += bsdfValue * mesh->getEmitter()->getRadiance() * fmaxf(0.0, Frame::cosTheta(bRec.wo)) / pdfOut;
            }else{
                pdfOut = NAN; //shadowRay被遮挡
            }
        }
        return result;
    }
    float getEmitterPdf(EmitterQueryRecord eRec, const Mesh* const emitterMesh) const
    {
        return emitterMesh->getEmitter()->pdf(eRec, emitterDPDF[0] / emitterMesh->getSurfaceArea());
    }
    EClassType getClassType() const { return EScene; }

private:
    std::vector<Mesh*> m_meshes;
    Integrator* m_integrator = nullptr;
    Sampler* m_sampler = nullptr;
    Camera* m_camera = nullptr;
    Accel* m_accel = nullptr;
    DiscretePDF emitterDPDF;
    std::vector<uint32_t> emitterMeshId {};
};

NORI_NAMESPACE_END
