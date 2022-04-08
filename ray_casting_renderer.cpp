/*
clang++ ray_casting_renderer.cpp -std=c++11 -o ray_casting_renderer && ./ray_casting_renderer
*/
#include<stdio.h>   // use right-hand coordinate
#include<math.h>    // time complexity: w*h*cnt_tri(for each T)*cnt_light(shade)*cnt_tri(shadow)
#include<vector>
#include<glm/glm.hpp>
#include<glm/ext.hpp>
#include<numeric>
using glm::vec3;
using std::vector;
using glm::u8;
class Ray {
public:
    Ray() {}
    Ray(const vec3& o, const vec3& d):org(o), dir(d) {}
    void print() { printf("org:(%.2f,%.2f,%.2f), dir:(%.2f,%.2f,%.2f)\n", org.x, org.y, org.z, dir.x, dir.y, dir.z); }
    vec3 org; vec3 dir;   // normalized
};
class BSDF {
public:
    BSDF():diff_color(vec3(0.0f)), spec_color(vec3(0.0f)), spec_gloss(0.0f) {}
    vec3 diff_color;
    vec3 spec_color;
    float spec_gloss;
};
class Triangle {
public:
    Triangle(const vec3& v1, const vec3& v2, const vec3& v3, const BSDF& _bsdf):bsdf(_bsdf) {
        vertices[0] = v1;   vertices[1] = v2;   vertices[2] = v3;
        vec3 normal = glm::normalize(glm::cross(vertices[2] - vertices[1], vertices[0] - vertices[2]));
        normals[0] = normals[1] = normals[2] = normal;
    } 
    vec3 vertices[3];
    vec3 normals[3];
    BSDF bsdf;
};
class Camera {
public:
    Camera(float fov, float zn, float a):deg_fov(fov), z_near(zn),aspect(a) {}
    float deg_fov;  // degree
    float z_near;
    float aspect;
};
class Light {
public:
    Light() {}
    Light(const vec3& o, const vec3 col):org(o), color(col) {}
    vec3 org;   vec3 color;
};
class Scene {
public:
    Scene() {}
    vector<Triangle> triangles;
    vector<Light> lights;
};
class Image {
public:
    Image(int _w, int _h):w(_w), h(_h) { colors = new glm::u8vec3[w*h]; }
    void fill_color(int x, int y, const glm::u8vec3& color) {
        int index = (h - 1 - y) * w + x;
        if(x >= 0 && y >= 0 && x < w && y < h)
            colors[index] = color;
    }
    void write(const char* filename) {
        FILE* f = fopen(filename, "w");
        fprintf(f, "P3\n%d %d\n255\n", w, h);
        for (int i = 0; i < w*h; i++) 
            fprintf(f, "%d %d %d ", colors[i].r, colors[i].g, colors[i].b);
        fclose(f);  
    }
    ~Image() { delete[] colors; }
    glm::u8vec3* colors;
    int w;  int h;
};
vec3 barycentric(const vec3& weights, const vec3 vecs[]) { return weights[0] * vecs[0] + weights[1] * vecs[1] + weights[2] * vecs[2]; }
int gamma_encode(float radiance) { return (int)(pow(std::min(1.0f, std::max(0.0f, radiance)), 1.0f/2.2f) * 255.0f); }
bool intersect(const Ray& r, const Triangle& tri, vec3& weights) { // pi: point of intersection, O + tD = b0P0 + b1P1 + b2P2
    vec3 O = r.org;             vec3 D = r.dir;
    vec3 P0 = tri.vertices[0];  vec3 P1 = tri.vertices[1];  vec3 P2 = tri.vertices[2];
    vec3 E1 = P1 - P0;
    vec3 E2 = P2 - P0;
    vec3 S = O - P0;
    vec3 S1 = glm::cross(D, E2);
    vec3 S2 = glm::cross(S, E1);
    float r_s1_dot_e1 = 1.0f / glm::dot(S1, E1);
    float t = r_s1_dot_e1 * glm::dot(S2, E2);
    float b1 = r_s1_dot_e1 * glm::dot(S1, S);
    float b2 = r_s1_dot_e1 * glm::dot(S2, D);
    float b0 = 1.0f - b1 - b2;
    if(t >= 0 && b1 >= 0 && b1 <= 1 && b2 >= 0 && b2 <= 1 && b0 >= 0 && b0 <= 1) {
        weights.x = b0; weights.y = b1; weights.z = b2;
        return true;
    }
    return false;
}
void cal_ray(Ray& r, int x, int y, int w, int h, const Camera& camera) {
    float half_height = camera.z_near * tanf(glm::radians(camera.deg_fov * 0.5f));
    float half_width = half_height * camera.aspect;
    r.org = vec3((((x + 0.5) / w) * 2.0f - 1.0f) * half_width, (((y + 0.5) / h) * 2.0f - 1.0f) * half_height, camera.z_near * -1);
    r.dir = glm::normalize(r.org);
}
bool visible(const Scene& scene, float dist2light, const vec3& p, const vec3& w_i) {
    float eps = 1e-4;
    Ray r(p + eps * w_i, w_i);
    dist2light -= eps;
    for(int i = 0; i < scene.triangles.size(); i++) {
        const Triangle& tri = scene.triangles[i];
        vec3 weights;
        if(intersect(r, tri, weights)) {
            vec3 pi = barycentric(weights, tri.vertices);
            if(glm::length(r.org - pi) < dist2light)return false;
        }
    }
    return true;
}
void shade(const Scene& scene, const Triangle& tri, 
           const vec3& pi, const vec3& normal, const vec3& w_o, vec3& radiance) {
    radiance = vec3(0.0f, 0.0f, 0.0f);
    for(int i = 0; i < scene.lights.size(); i++) {
        const Light& light = scene.lights[i];
        vec3 w_i = light.org - pi;
        float dist2light = glm::length(w_i);
        w_i = w_i / dist2light;
        if(visible(scene, dist2light, pi, w_i)) {
            float recip_square_r = 1.0f / (dist2light*dist2light);
            // diffuse: k_diff * c_diff * c_light * 1/r^2 * max(0, dot(n, w_i))
            float k_diff = 5.0f;
            radiance += k_diff * tri.bsdf.diff_color * light.color * recip_square_r
                        * std::max(0.0f, glm::dot(normal, w_i));
            // specular: k_spec * c_spec * c_light * 1/r^2 * max(0,dot(h, n))^spec_gloss
            float k_spec = 2.0f;
            vec3 h = glm::normalize(w_o + w_i);
            radiance += k_spec * tri.bsdf.spec_color * light.color * recip_square_r
                        * pow(std::max(0.0f, glm::dot(h, normal)), tri.bsdf.spec_gloss);
        }
    }
}
void ray_trace(Image& image, const Scene& scene, const Camera& camera) {
    int w = image.w, h = image.h;
    for(int y = 0; y < h; y++) {
        for(int x = 0; x < w; x++) {
            Ray r;
            cal_ray(r, x, y, w, h, camera);
            float closest = std::numeric_limits<float>::infinity();
            vec3 radiance;
            for(int i = 0; i < scene.triangles.size(); i++) {
                const Triangle& tri = scene.triangles[i];
                vec3 weights;
                if(intersect(r, tri, weights)) {
                    vec3 pi = barycentric(weights, tri.vertices);
                    vec3 normal = barycentric(weights, tri.normals);
                    float dist = glm::distance(pi, r.org);
                    if(dist < closest) {
                        closest = dist;
                        shade(scene, tri, pi, normal, -1.0f * r.dir, radiance);
                        glm::u8vec3 col(gamma_encode(radiance.x), gamma_encode(radiance.y), gamma_encode(radiance.z));
                        image.fill_color(x, y, col);
                    }
                }
            }
        }
    }
    image.write("ray_casting_renderer.ppm");
}
void init_scene(Scene& scene) {
    BSDF bsdf;
    bsdf.diff_color = vec3(0.05f, 0.4f, 0.1f);
    bsdf.spec_color = vec3(1.0f, 1.0f, 1.0f);
    bsdf.spec_gloss = 64.0f;

    Triangle tri(vec3(  -0.8f + 0.2f,   2.0f, -4.0f), 
                 vec3(  -1.6f + 0.2f,  -1.0f, -2.0f), 
                 vec3(   1.5f + 0.2f,  -1.0f, -4.0f), 
                 bsdf);
    scene.triangles.push_back(tri);

    BSDF mat_plane;
    mat_plane.diff_color = 4.0f * vec3(1.0f);
    mat_plane.spec_color = vec3(0.0f, 0.0f, 0.0f);
    mat_plane.spec_gloss = 0.0f;

    float half_a = 30.0f, bottom_y = -3.0f, top_y = 50.0f;
    vec3 p0 = vec3( half_a,   bottom_y,   -half_a);
    vec3 p1 = vec3(-half_a,   bottom_y,    half_a);
    vec3 p2 = vec3( half_a,   bottom_y,    half_a);
    vec3 p3 = vec3(-half_a,   bottom_y,   -half_a);
    Triangle plane1(p0, p1, p2, mat_plane);
    Triangle plane2(p1, p0, p3, mat_plane);
    scene.triangles.push_back(plane1);
    scene.triangles.push_back(plane2);

    p0 = vec3(-half_a,  bottom_y,   -half_a);
    p1 = vec3( half_a,  bottom_y,   -half_a);
    p2 = vec3( half_a,     top_y,   -half_a);
    p3 = vec3(-half_a,     top_y,   -half_a);
    plane1 = Triangle(p0, p1, p2, mat_plane);
    plane2 = Triangle(p0, p2, p3, mat_plane);
    scene.triangles.push_back(plane1);
    scene.triangles.push_back(plane2);

    Light light(vec3(0.0f, 1.0f, -2.0f), vec3(1.0f, 1.0f, 1.0f));
    scene.lights.push_back(light);
    return;
}
int main() {
    int w = 1024, h = 1024;
    Image image(w, h);
    Camera camera(90.0f, 1.0f, (float)w / h);
    Scene scene;
    init_scene(scene);
    ray_trace(image, scene, camera);
    return 0;
}