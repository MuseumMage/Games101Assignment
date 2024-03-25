//
// Created by Göksu Güvendiren on 2019-05-14.
//

#include "Scene.hpp"


void Scene::buildBVH() {
    printf(" - Generating BVH...\n\n");
    this->bvh = new BVHAccel(objects, 1, BVHAccel::SplitMethod::NAIVE);
}

Intersection Scene::intersect(const Ray &ray) const
{
    return this->bvh->Intersect(ray);
}

void Scene::sampleLight(Intersection &pos, float &pdf) const
{
    float emit_area_sum = 0;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        if (objects[k]->hasEmit()){
            emit_area_sum += objects[k]->getArea();
        }
    }
    float p = get_random_float() * emit_area_sum;
    emit_area_sum = 0;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        if (objects[k]->hasEmit()){
            emit_area_sum += objects[k]->getArea();
            if (p <= emit_area_sum){
                objects[k]->Sample(pos, pdf);
                break;
            }
        }
    }
}

bool Scene::trace(
        const Ray &ray,
        const std::vector<Object*> &objects,
        float &tNear, uint32_t &index, Object **hitObject)
{
    *hitObject = nullptr;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        float tNearK = kInfinity;
        uint32_t indexK;
        Vector2f uvK;
        if (objects[k]->intersect(ray, tNearK, indexK) && tNearK < tNear) {
            *hitObject = objects[k];
            tNear = tNearK;
            index = indexK;
        }
    }
    return (*hitObject != nullptr);
}

// Implementation of Path Tracing
Vector3f Scene::castRay(const Ray &ray, int depth) const
{
    // 可参考 https://zhuanlan.zhihu.com/p/488882096
    // Implement Path Tracing Algorithm here
    // init intersection and w0
    auto intersection = intersect(ray);
    if (!intersection.happened) {
        return {};
    }

    // Emission
    if (intersection.m->hasEmission()) {
        return intersection.m->getEmission();
    }

    // init
    // w0的方向好像不影响？
    // 原因应该是eval和pdf里面第一个参数没有被用到
    auto w0 = -ray.direction;
    auto L_dir = Vector3f();
    auto L_indir = Vector3f();

    // 1. from light source
    // Uniformly sample the light at x` (pdf_light = 1 / A)
    // L_dir = L_i * f_r * cos θ * cos θ` / |x` - intersection|^2 / pdf_light
    float pdf_light = 0;
    Intersection hit_light;
    sampleLight(hit_light, pdf_light);
    auto p = intersection.coords;
    auto x = hit_light.coords;
    auto ws_unnorm = x - p;
    auto ws = ws_unnorm.normalized();
    auto nn = hit_light.normal;

    // Shoot a ray from intersection to x
    Ray block_ray(p, ws);
    // Check if the ray is blocked
    Intersection block_intersect = intersect(block_ray);
//    if (!block_intersect.happened)
    if (block_intersect.distance - ws_unnorm.norm() > -0.005)
    {
        auto L_i = hit_light.emit;
        auto f_r = intersection.m->eval(w0, ws, intersection.normal);
        auto cos_theta = std::max(0.0f, dotProduct(intersection.normal, ws));
        auto cos_theta_prime = std::max(0.0f, dotProduct(nn, -ws));
        auto r2 = dotProduct(ws_unnorm, ws_unnorm);
        L_dir = L_i * f_r * cos_theta * cos_theta_prime / r2 / pdf_light;
    }

    // 2. from indirect light
    //  L_indir = shade(q, wi) * f_r * cos_theta / pdf_hemi / P_RR
    // RussianRoulette Test
    float ksi = get_random_float();
    if (ksi > RussianRoulette)
    {
        return L_dir + L_indir;
    }

    auto wi = (intersection.m->sample(w0, intersection.normal)).normalized();
    auto secondary_ray = Ray(p, wi);
    auto secondary_inter= intersect(secondary_ray);
    if (secondary_inter.happened && !secondary_inter.m->hasEmission())
    {
        auto f_r = intersection.m->eval(w0, wi, intersection.normal);
        auto cos_theta = std::max(0.0f, dotProduct(intersection.normal, wi));
        auto pdf_hemi = intersection.m->pdf(w0, wi, intersection.normal);
        L_indir = castRay(secondary_ray, depth + 1) * f_r * cos_theta / pdf_hemi / RussianRoulette;
    }


    return L_dir + L_indir;
}