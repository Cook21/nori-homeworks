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

#include "nori/color.h"
#include "nori/common.h"
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

    Color3f sampleEmitter(Sampler* sampler, const Point3f& shadingPoint, Vector3f& outDir) const
    {
        if (emitterSurfaceAreaDPDF.size() > 0) {
            float sampleFacePdf;
            auto id = emitterSurfaceAreaDPDF.sample(sampler->next1D(), sampleFacePdf);
            Mesh* mesh = m_meshes[emitterIdx[id]];
            Point3f lightSamplePos;
            Vector3f lightSamplePosNormal;
            float SamplePosPdf;
            mesh->sample(sampler, lightSamplePos, lightSamplePosNormal, SamplePosPdf);
            outDir = lightSamplePos - shadingPoint;
            float distanceSquared = outDir.dot(outDir);
            //算完距离之后单位化
            outDir.normalize();
            auto secondaryRay = Ray3f(shadingPoint, outDir);
            Intersection shadowRayIts;
            this->shadowrayIntersect(secondaryRay, shadowRayIts);
            if (shadowRayIts.t * shadowRayIts.t >= distanceSquared - Epsilon) {
                return mesh->getEmitter()->sample(-outDir, lightSamplePosNormal, distanceSquared) / (SamplePosPdf * sampleFacePdf);
            } else {
                return Color3f(0.0f);
            }
        } else {
            return Color3f(0.f);
        }
    }

    EClassType getClassType() const { return EScene; }

private:
    std::vector<Mesh*> m_meshes;
    Integrator* m_integrator = nullptr;
    Sampler* m_sampler = nullptr;
    Camera* m_camera = nullptr;
    Accel* m_accel = nullptr;
    DiscretePDF emitterSurfaceAreaDPDF;
    std::vector<uint32_t> emitterIdx {};
};

NORI_NAMESPACE_END
