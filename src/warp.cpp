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

#include <cmath>
#include <nori/frame.h>
#include <nori/vector.h>
#include <nori/warp.h>

NORI_NAMESPACE_BEGIN

Point2f Warp::squareToUniformSquare(const Point2f& sample)
{
    return sample;
}

float Warp::squareToUniformSquarePdf(const Point2f& sample)
{
    return ((sample.array() >= 0).all() && (sample.array() <= 1).all()) ? 1.0f : 0.0f;
}

Point2f Warp::squareToTent(const Point2f& sample)
{
    auto inverseCDF = [](float x) -> float {
        const float sqrtOf2 = sqrt(2);
        if (0 <= x && x < .5) {
            return sqrtOf2 * sqrt(x) - 1.0;
        } else if (.5 <= x && x <= 1.) {
            return -sqrtOf2 * sqrt(1.0 - x) + 1.0;
        } else {
            return 0.;
        }
    };
    return Point2f(inverseCDF(sample[0]), inverseCDF(sample[1]));
}

float Warp::squareToTentPdf(const Point2f& p)
{
    auto PDF = [](float x) -> float {
        if (-1. <= x && x < 0.) {
            return x + 1.0;
        } else if (0. <= x && x <= 1.) {
            return 1.0 - x;
        } else {
            return 0.;
        }
    };
    return PDF(p[0]) * PDF(p[1]);
}

Point2f Warp::squareToUniformDisk(const Point2f& sample)
{
    auto uniformSampleDisk = [&sample]() -> Point2f {
        float r = sqrt(sample[0]);
        float theta = 2.0 * M_PI * sample[1];
        return Point2f(r * cos(theta), r * sin(theta));
    };
    auto concentricSampleDisk = [&sample]() -> Point2f {
        const float pi_over_4 = M_PI / 4;
        const float pi_over_2 = M_PI / 2;
        //map sample from [0,1] to [-1,1]
        Point2f newSample = 2.f * sample - Vector2f(1, 1);
        if (newSample[0] == 0 && newSample[1] == 0) {
            return Point2f(0, 0);
        } else {
            float r, theta;
            if (abs(newSample[0]) > abs(newSample[1])) {
                r = newSample[0];
                theta = pi_over_4 * newSample[1] / newSample[0];
            } else {
                r = newSample[1];
                theta = pi_over_2 - (pi_over_4 * newSample[0] / newSample[1]);
            }
            return r * Point2f(cos(theta), sin(theta));
        }
    };

    return concentricSampleDisk();
}

float Warp::squareToUniformDiskPdf(const Point2f& p)
{
    float square_of_r = p[0] * p[0] + p[1] * p[1];
    if (square_of_r <= 1) {
        return INV_PI;
    } else {
        return 0.0;
    }
}

Vector3f Warp::squareToUniformSphere(const Point2f& sample)
{
    float temp1 = 2. * sqrt((1.0 - sample[0]) * sample[0]);
    float temp2 = 2. * M_PI * sample[1];
    float x = cos(temp2) * temp1;
    float y = sin(temp2) * temp1;
    float z = 1. - 2. * sample[0];

    return Vector3f(x, y, z);
}

float Warp::squareToUniformSpherePdf(const Vector3f& v)
{
    return INV_FOURPI;
}

Vector3f Warp::squareToUniformHemisphere(const Point2f& sample)
{
    float z = sample[0];
    float temp1 = sqrt(1.0 - z * z);
    float temp2 = 2. * M_PI * sample[1];
    float x = cos(temp2) * temp1;
    float y = sin(temp2) * temp1;

    return Vector3f(x, y, z);
}

float Warp::squareToUniformHemispherePdf(const Vector3f& v)
{
    if (v[2] >= 0) {
        return INV_TWOPI;
    } else {
        return 0.0;
    }
}

Vector3f Warp::squareToCosineHemisphere(const Point2f& sample)
{
    Point2f projectionPoint = squareToUniformDisk(sample);
    float x = projectionPoint[0];
    float y = projectionPoint[1];
    float z = sqrt(1. - x * x - y * y);
    return Vector3f(x, y, z);
}

float Warp::squareToCosineHemispherePdf(const Vector3f& v)
{
    if (v[2] >= 0) {
        float r = sqrt(v[0] * v[0] + v[1] * v[1]);
        float phi = atan(r / v[2]);
        return cos(phi) * INV_PI;
    } else {
        return 0;
    }
}

Vector3f Warp::squareToBeckmann(const Point2f& sample, float alpha)
{
    float theta = 2. * M_PI * sample[0];
    float tanPhi2 = -alpha*alpha * log(1 - sample[1]);
    float cosPhi = sqrt(1 / (1 + tanPhi2));
    float sinPhi = sqrt(1. - cosPhi * cosPhi);
    return Vector3f(cos(theta) * sinPhi, sin(theta) * sinPhi, cosPhi);
}

float Warp::squareToBeckmannPdf(const Vector3f& m, float alpha)
{
    if (m[2] >= 0) {
        float r = sqrt(m[0] * m[0] + m[1] * m[1]);
        float phi = atan(r / m[2]);
        float squareOfAlpha = alpha * alpha;
        return INV_PI * exp(-tan(phi) * tan(phi) / squareOfAlpha) / squareOfAlpha / pow(cos(phi), 3);
    } else {
        return 0;
    }
}

NORI_NAMESPACE_END
