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
#include "nori/mesh.h"
#include <Eigen/Geometry>
#include <chrono>
#include <cstdint>
#include <functional>
#include <new>
#include <nori/accel.h>

NORI_NAMESPACE_BEGIN

void Accel::addMesh(Mesh* mesh)
{
    if (m_mesh)
        throw NoriException("Accel: only a single mesh is supported!");
    m_mesh = mesh;
    m_bbox = m_mesh->getBoundingBox();
}

void Accel::build()
{
    std::function<OctTreeNode*(TBoundingBox<Point3f>*, std::vector<uint32_t>*, nori::Mesh*, int depth)> buildTree;
    buildTree = [&buildTree](TBoundingBox<Point3f>* box, std::vector<uint32_t>* triangles, nori::Mesh* mesh, int depth) -> OctTreeNode* {
        const int maxDepth = 10;
        const int maxTriangleCount = 16;
        if (triangles->size() == 0 || !box->hasVolume()) {
            return nullptr;
        } else if (triangles->size() <= maxTriangleCount || depth > maxDepth) {
            return new OctTreeNode(box, triangles);
        } else {
            auto* parent = new OctTreeNode(box, nullptr);
            auto center = box->getCenter();
            std::array<TBoundingBox<Point3f>*, 8> childBBox;
            std::array<std::vector<uint32_t>*, 8> childTriangles;
            for (int i = 0; i < 8; i++) {
                auto corner = box->getCorner(i);
                childBBox[i] = new TBoundingBox<Point3f>(center);
                childBBox[i]->expandBy(corner);
                childTriangles[i] = new std::vector<uint32_t>;
            }
            for (uint32_t triangle : *triangles) {
                TBoundingBox<Point3f> triangleBBox = mesh->getBoundingBox(triangle);
                for (int i = 0; i < 8; i++) {
                    if (childBBox[i]->overlaps(triangleBBox)) {
                        childTriangles[i]->emplace_back(triangle);
                    }
                }
            }
            for (int i = 0; i < 8; i++) {
                parent->child[i] = buildTree(childBBox[i], childTriangles[i], mesh, depth + 1);
            }
            delete triangles;

            return parent;
        }
    };
    using namespace std::chrono;
    auto start = high_resolution_clock::now();
    //put all triangles into list
    auto* triangles = new std::vector<uint32_t>(m_mesh->getTriangleCount());
    for (uint32_t idx = 0; idx < m_mesh->getTriangleCount(); ++idx) {
        (*triangles)[idx] = idx;
    }
    m_octTree = buildTree(new TBoundingBox<Point3f>(m_bbox), triangles, m_mesh, 0);
    auto end = high_resolution_clock::now();
    std::cout << "OctTree build time:" << duration_cast<milliseconds>(end - start).count() << "ms\n";
}
void Accel::sortTree(const Point3f& origin)
{
    std::function<void(OctTreeNode & node, const Point3f& origin)> sortOctTree;
    sortOctTree = [&sortOctTree](OctTreeNode& node, const Point3f origin) -> void {
        std::sort(node.child.begin(), node.child.end(), [&origin](OctTreeNode* a, OctTreeNode* b) -> bool {
            if (a == b) {
                return false;
            } else if (a == nullptr) {
                return false;
            } else if (b == nullptr) {
                return true;
            } else {
                return a->boundingBox->squaredDistanceTo(origin) < b->boundingBox->squaredDistanceTo(origin);
            }
        });
        for (auto* i : node.child) {
            if (i != nullptr) {
                sortOctTree(*i, origin);
            }
        }
    };
    sortOctTree(*m_octTree, origin);
}

bool Accel::rayIntersect(const Ray3f& ray_, Intersection& its, bool shadowRay) const
{
    auto bruteForceSearch = [](Ray3f& ray, Intersection& its, bool shadowRay, nori::Mesh* mesh) -> std::pair<bool, uint32_t> {
        bool didFound = false;
        uint32_t itsTriangle = (uint32_t)-1;
        for (uint32_t idx = 0; idx < mesh->getTriangleCount(); ++idx) {
            float u, v, t;
            if (mesh->rayIntersect(idx, ray, u, v, t)) {
                /* An intersection was found! Can terminate
               immediately if this is a shadow ray query */
                if (shadowRay) {
                    return { true, itsTriangle };
                } else {
                    ray.maxt = its.t = t;
                    its.uv = Point2f(u, v);
                    its.mesh = mesh;
                    didFound = true;
                    itsTriangle = idx;
                }
            }
        }
        return { didFound, itsTriangle };
    };
    std::function<std::pair<bool, uint32_t>(const OctTreeNode& node, Ray3f& ray, Intersection& its, nori::Mesh* mesh)> octTreeSearch;
    octTreeSearch = [&octTreeSearch](const OctTreeNode& node, Ray3f& ray, Intersection& its, nori::Mesh* mesh) -> std::pair<bool, uint32_t> {
        bool didFound = false;
        uint32_t itsTriangle = (uint32_t)-1;
        if (node.boundingBox->rayIntersect(ray)) {

            if (node.isLeaf()) {
                for (auto triangle : *(node.triangles)) {
                    float u, v, t;
                    if (mesh->rayIntersect(triangle, ray, u, v, t)) {
                        ray.maxt = its.t = t;
                        its.uv = Point2f(u, v);
                        its.mesh = mesh;
                        didFound = true;
                        itsTriangle = triangle;
                    }
                }
            } else {
                for (int i = 0; i < 8; i++) {
                    if (node.child[i] != nullptr) {
                        auto result = octTreeSearch(*(node.child[i]), ray, its, mesh);
                        if (result.first) {
                            didFound = result.first;
                            itsTriangle = result.second;
                            break; //Only when the octTree nodes have been sort can we break
                        }
                    }
                }
            }
        }
        return { didFound, itsTriangle };
    };

    std::function<std::pair<bool, uint32_t>(const OctTreeNode& node, const Ray3f& ray, nori::Mesh* mesh)> octTreeSearchForShadowRay;
    octTreeSearchForShadowRay = [&octTreeSearchForShadowRay](const OctTreeNode& node, const Ray3f& ray, nori::Mesh* mesh) -> std::pair<bool, uint32_t> {
        bool didFound = false;
        uint32_t itsTriangle = (uint32_t)-1;
        if (node.boundingBox->rayIntersect(ray)) {
            if (node.isLeaf()) {
                for (auto triangle : *(node.triangles)) {
                    float u, v, t;
                    if (mesh->rayIntersect(triangle, ray, u, v, t)) {
                        didFound = true;
                        break;
                    }
                }
            } else {
                for (int i = 0; i < 8; i++) {
                    if (node.child[i] != nullptr) {
                        auto result = octTreeSearchForShadowRay(*(node.child[i]), ray, mesh);
                        if (result.first) {
                            didFound = true;
                            break;
                        }
                    }
                }
            }
        }
        return { didFound, itsTriangle };
    };

    Ray3f ray(ray_); /// Make a copy of the ray (we will need to update its '.maxt' value), maxt acts as a z-buffer

    //1. original method: Brute force search through all triangles
    //auto [didFoundIntersection, itsTriangle] = bruteForceSearch(ray, its, shadowRay, m_mesh);

    //2. acclerated method

    bool didFoundIntersection = false;
    uint32_t itsTriangle = (uint32_t)-1;
    if (shadowRay) {
        auto result = octTreeSearchForShadowRay(*m_octTree, ray, m_mesh);
        didFoundIntersection = result.first;
    } else {
        auto result = octTreeSearch(*m_octTree, ray, its, m_mesh);
        didFoundIntersection = result.first;
        itsTriangle = result.second;
    }

    if (shadowRay) {
        return didFoundIntersection;
    } else {
        if (didFoundIntersection) {
            /* At this point, we now know that there is an intersection,
           and we know the triangle index of the closest such intersection.

           The following computes a number of additional properties which
           characterize the intersection (normals, texture coordinates, etc..)
        */

            /* Find the barycentric coordinates */
            Vector3f bary;
            bary << 1 - its.uv.sum(), its.uv;

            /* References to all relevant mesh buffers */
            const Mesh* mesh = its.mesh;
            const MatrixXf& V = mesh->getVertexPositions();
            const MatrixXf& N = mesh->getVertexNormals();
            const MatrixXf& UV = mesh->getVertexTexCoords();
            const MatrixXu& F = mesh->getIndices();

            /* Vertex indices of the triangle */
            uint32_t idx0 = F(0, itsTriangle), idx1 = F(1, itsTriangle), idx2 = F(2, itsTriangle);

            Point3f p0 = V.col(idx0), p1 = V.col(idx1), p2 = V.col(idx2);

            /* Compute the intersection positon accurately
           using barycentric coordinates */
            its.p = bary.x() * p0 + bary.y() * p1 + bary.z() * p2;

            /* Compute proper texture coordinates if provided by the mesh */
            if (UV.size() > 0)
                its.uv = bary.x() * UV.col(idx0) + bary.y() * UV.col(idx1) + bary.z() * UV.col(idx2);

            /* Compute the geometry frame */
            its.geoFrame = Frame((p1 - p0).cross(p2 - p0).normalized());

            if (N.size() > 0) {
                /* Compute the shading frame. Note that for simplicity,
               the current implementation doesn't attempt to provide
               tangents that are continuous across the surface. That
               means that this code will need to be modified to be able
               use anisotropic BRDFs, which need tangent continuity */

                its.shFrame = Frame(
                    (bary.x() * N.col(idx0) + bary.y() * N.col(idx1) + bary.z() * N.col(idx2)).normalized());
            } else {
                its.shFrame = its.geoFrame;
            }
        }

        return didFoundIntersection;
    }
}

NORI_NAMESPACE_END
