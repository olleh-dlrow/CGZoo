/*
c++ rasterization_renderer.cpp -std=c++11 -o rasterization_renderer && ./rasterization_renderer
*/
#include<stdio.h>
#include<vector>
#include<numeric>
#include<glm/ext.hpp>
#include<glm/glm.hpp>
using std::vector;
using glm::vec2;
using glm::vec3;
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
    Image(int _w, int _h):w(_w), h(_h) { 
        colors = new glm::u8vec3[w*h];
        for(int i = 0; i < w*h; i++) {
            colors[i] = glm::u8vec3(255, 255, 255);
        }
        dist_buffer = new float[w*h]; 
        for(int i = 0; i < w*h; i++) {
            dist_buffer[i] = std::numeric_limits<float>::infinity();
        }
    }
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
    ~Image() { delete[] colors; delete[] dist_buffer; }
    glm::u8vec3* colors;
    float* dist_buffer;
    int w;  int h;
};
float line_dist2D(const vec2& p, const vec2& a, const vec2& b) {
    vec2 n = vec2(a.y - b.y, b.x - a.x);    // xAB
    float d = a.x * b.y - a.y * b.x;
    return (glm::dot(n, p) + d) / glm::length(n);
}
float cal_bary2D(const vec2& a, const vec2& b, const vec2& c, const vec2& p) { return line_dist2D(p, b, c) / line_dist2D(a, b, c); }
vec2 barycentric2D(const vec2& weights, const vec2 vecs[]) { return weights[0] * vecs[0] + weights[1] + vecs[1]; }
vec3 barycentric3D(const vec3& weights, const vec3 vecs[]) { return weights[0] * vecs[0] + weights[1] * vecs[1] + weights[2] * vecs[2]; }
float barycentric3D(const vec3& weights, const float fs[]) {return weights[0] * fs[0] + weights[1] * fs[1] + weights[2] * fs[2]; }
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
            vec3 pi = barycentric3D(weights, tri.vertices);
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
            // diffuse
            // k_diff * c_diff * c_light * 1/r^2 * max(0, dot(n, w_i))
            float k_diff = 5.0f;
            radiance += k_diff * tri.bsdf.diff_color * light.color * 1.0f / (dist2light*dist2light)
                        * std::max(0.0f, glm::dot(normal, w_i));
            
            float k_spec = 2.0f;
            vec3 h = glm::normalize(w_o + w_i);
            radiance += k_spec * tri.bsdf.spec_color * light.color * 1.0f / (dist2light*dist2light)
                        * pow(std::max(0.0f, glm::dot(h, normal)), tri.bsdf.spec_gloss);
        }
    }
}
// 只有三角形的z大于zNear时才能调用，否则计算出的投影点是错误的！！！
void bounding_box(const Camera& camera, const Triangle& tri, glm::i32vec2 proj_vecs[],
    int w, int h, int& x_min, int& y_min, int& x_max, int& y_max) {
    x_min = w - 1;  x_max = 0;
    y_min = h - 1;  y_max = 0;
    for(int i = 0; i < 3; i++) {
        const vec3& v = tri.vertices[i];
        vec3 v_per = v * -1.0f / v.z;
        float half_h = glm::tan(glm::radians(camera.deg_fov * 0.5f));
        float half_w = half_h * camera.aspect;
        float y = (v_per.y / half_h + 1) * 0.5f * (h - 1);
        float x = (v_per.x / half_w + 1) * 0.5f * (w - 1);
        proj_vecs[i] = glm::i32vec2(x, y);
        x_min = glm::min(int(x - 0.5f), x_min);
        x_max = glm::max(int(x + 0.5f), x_max);
        y_min = glm::min(int(y - 0.5f), y_min);
        y_max = glm::max(int(y + 0.5f), y_max);
    }
    x_min = glm::max(    0, x_min);
    x_max = glm::min(w - 1, x_max);
    y_min = glm::max(    0, y_min);
    y_max = glm::min(h - 1, y_max);    
}
void rasterization(Image& image, const Scene& scene, const Camera& camera) {
    int w = image.w, h = image.h;
    for(int i = 0; i < scene.triangles.size(); i++) {
        const Triangle& tri = scene.triangles[i];
        int x_min = 0, y_min = 0, x_max = w - 1, y_max = h - 1;
        glm::i32vec2 proj_vecs[3];
        bounding_box(camera, tri, proj_vecs, w, h, x_min, y_min, x_max, y_max);
        vec3 vertexPw[3];
        vec3 vertexNw[3];
        float vertexW[3];
        for(int j = 0; j < 3; j++) {
            vertexW [j] = -1.0f / tri.vertices[j].z;
            vertexPw[j] = tri.vertices[j] * vertexW[j];
            vertexNw[j] = tri.normals [j] * vertexW[j];
        }
        // printf("x:%d,%d y:%d,%d\n", x_min, x_max, y_min, y_max);
        for(int y = y_min; y <= y_max; y++) {
            for(int x = x_min; x <= x_max; x++) {
                Ray r;
                cal_ray(r, x, y, w, h, camera);
                vec3 weights = vec3(
                    cal_bary2D(proj_vecs[0], proj_vecs[1], proj_vecs[2], vec2(x, y)),
                    cal_bary2D(proj_vecs[1], proj_vecs[2], proj_vecs[0], vec2(x, y)),
                    cal_bary2D(proj_vecs[2], proj_vecs[0], proj_vecs[1], vec2(x, y))
                );
                if(weights[0] >= 0 && weights[0] <= 1 &&
                   weights[1] >= 0 && weights[1] <= 1 &&
                   weights[2] >= 0 && weights[2] <= 1   ) {
                    vec3 inter_p      = barycentric3D(weights, vertexPw) / barycentric3D(weights, vertexW);
                    vec3 inter_normal = barycentric3D(weights, vertexNw) / barycentric3D(weights, vertexW);
                    inter_normal = glm::normalize(inter_normal);
                    float dist = glm::abs(inter_p.z);
                    if(dist < image.dist_buffer[(h - 1 - y)*w + x]) {
                        image.dist_buffer[(h - 1 - y)*w + x] = dist;
                        vec3 radiance;
                        shade(scene, tri, inter_p, inter_normal, -1.0f*r.dir, radiance);
                        glm::u8vec3 col(gamma_encode(radiance.r), gamma_encode(radiance.g), gamma_encode(radiance.b));
                        image.fill_color(x, y, col);
                    }
                }
            }
        }
    }
    image.write("rasterization_renderer.ppm");
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
    vec3 p1 = vec3(-half_a,   bottom_y,     -1.1f);
    vec3 p2 = vec3( half_a,   bottom_y,     -1.1f);
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
    rasterization(image, scene, camera);
    return 0;
}