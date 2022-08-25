/*
c++ rendering_pipeline.cpp -std=c++17 -o rendering_pipeline && ./rendering_pipeline
*/
#include<stdio.h>
#include<vector>
#include<numeric>
#include<functional>
#include<glm/ext.hpp>
#include<glm/glm.hpp>
#include"thread_pool.h"
using std::vector;
using std::function;
using glm::vec2;
using glm::vec3;
using glm::vec4;
using glm::u8vec3;
using glm::mat4;
class Camera {
public:
    Camera(float fov, float zn, float zf, float a):deg_fov(fov), z_near(zn), z_far(zf), aspect(a) {}
    float deg_fov;  float z_near;   float z_far;    float aspect;
};
class BSDF {
public:
    BSDF():diff_color(vec3(0.0f)), spec_color(vec3(0.0f)), spec_gloss(0.0f) {}
    vec3 diff_color;
    vec3 spec_color;
    float spec_gloss;
};
class TriangleAttribute {
public:
    vec4 vertices[3];  // store homogeneous coordinates and ndc coordinates
    vec3 normals [3];  // 
    vec3 wd_pos  [3];  // world pos
    BSDF bsdf;
};
class Fragment {
public:  
    Fragment() {}
    vec3  wd_pos;
    vec3  normal;
    BSDF  bsdf;
};
class Triangle {
public:
    Triangle() {}
    Triangle(const vec3& v1, const vec3& v2, const vec3& v3, const BSDF& _bsdf):bsdf(_bsdf) {
        vertices[0] = v1;   vertices[1] = v2;   vertices[2] = v3;
        cal_normal();
    }
    void cal_normal() {
        vec3 normal = glm::normalize(glm::cross(vertices[2] - vertices[1], vertices[0] - vertices[2]));
        normals[0] = normals[1] = normals[2] = normal;
    }
    vec3 vertices[3];
    vec3 normals[3];
    BSDF bsdf;
};
class Polygon {
public:
    Polygon() {}
    Polygon(const Triangle& tri) {
        vertices.push_back(vec4(tri.vertices[0], 1));
        vertices.push_back(vec4(tri.vertices[1], 1));
        vertices.push_back(vec4(tri.vertices[2], 1));
    }
    void apply(const mat4& M) { for(int i = 0; i < vertices.size(); i++)vertices[i] = M * vertices[i]; }
    vector<vec4> vertices;
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
        z_buffer = new float[w*h]; 
        for(int i = 0; i < w*h; i++) {
            z_buffer[i] = std::numeric_limits<float>::infinity();
        }
    }
    inline int get_index(int x, int y) const { return (h - 1 - y) * w + x; }
    void get_xy(int index, int& x, int& y) const { x = index%w; y = h - 1 - index/w; }
    void fill_color(int x, int y, const glm::u8vec3& color) {
        int index = get_index(x, y);
        if(x >= 0 && y >= 0 && x < w && y < h) colors[index] = color;
    }
    void write(const char* filename) {
        FILE* f = fopen(filename, "w");
        fprintf(f, "P3\n%d %d\n255\n", w, h);
        for (int i = 0; i < w*h; i++) 
            fprintf(f, "%d %d %d ", colors[i].r, colors[i].g, colors[i].b);
        fclose(f);  
    }
    ~Image() { delete[] colors; delete[] z_buffer; }
    glm::u8vec3* colors;
    float* z_buffer;
    int w;  int h;
};
inline bool in_same_side(float s1, float s2){ return s1 * s2 >= 0; }
Polygon sutherland_hodgman(const Polygon& input, 
                           function<float(const vec4& v)> sign, 
                           const vec4& inner,
                           function<float(const vec4& s, const vec4& e)> cal_t) {
    Polygon res;
    for(int i = 0; i < input.vertices.size(); i++) {
        const vec4& s = input.vertices[i];
        const vec4& e = input.vertices[(i + 1) % input.vertices.size()];
        float inner_sign = sign(inner);
        float s_sign     = sign(s);
        float e_sign     = sign(e);
        if(!in_same_side(s_sign, inner_sign) && !in_same_side(e_sign, inner_sign)) {// 1. s out, e out: output None
            continue;
        } else if(in_same_side(s_sign, inner_sign) && in_same_side(e_sign, inner_sign)) {// 2. s in, e in: output e
            res.vertices.push_back(e);
        } else if(!in_same_side(s_sign, inner_sign) && in_same_side(e_sign, inner_sign)) {// 3. s out, e in: output i, e
            float t = cal_t(s, e);
            vec4 intersect = (1 - t) * s + t * e;
            res.vertices.push_back(intersect);
            res.vertices.push_back(e);
        } else {// 4. s in, e out: output i
            float t = cal_t(s, e);
            vec4 intersect = (1 - t) * s + t * e;
            res.vertices.push_back(intersect);
        }
    }
    return res;
}
Polygon homo_clip(const Polygon& input) {
    Polygon poly = input;
    vec4 inner = vec4(0, 0, 0, 1);
    // w + x = 0
    poly = sutherland_hodgman(poly, 
        [=](const vec4& v)->float{ return v.w + v.x; }, inner,
        [=](const vec4& s, const vec4& e)->float{ return ( s.w + s.x) / (( s.w + s.x) - ( e.w + e.x)); });
    // w - x = 0
    poly = sutherland_hodgman(poly, 
        [=](const vec4& v)->float{ return v.w - v.x; }, inner,
        [=](const vec4& s, const vec4& e)->float{ return (-s.w + s.x) / ((-s.w + s.x) - (-e.w + e.x)); });
    // w + y = 0
    poly = sutherland_hodgman(poly, 
        [=](const vec4& v)->float{ return v.w + v.y; }, inner,
        [=](const vec4& s, const vec4& e)->float{ return ( s.w + s.y) / (( s.w + s.y) - ( e.w + e.y)); });
    // w - y = 0
    poly = sutherland_hodgman(poly, 
        [=](const vec4& v)->float{ return v.w - v.y; }, inner,
        [=](const vec4& s, const vec4& e)->float{ return (-s.w + s.y) / ((-s.w + s.y) - (-e.w + e.y)); });
    // w + z = 0
    poly = sutherland_hodgman(poly, 
        [=](const vec4& v)->float{ return v.w + v.z; }, inner,
        [=](const vec4& s, const vec4& e)->float{ return ( s.w + s.z) / (( s.w + s.z) - ( e.w + e.z)); });
    // w - z = 0
    poly = sutherland_hodgman(poly, 
        [=](const vec4& v)->float{ return v.w - v.z; }, inner,
        [=](const vec4& s, const vec4& e)->float{ return (-s.w + s.z) / ((-s.w + s.z) - (-e.w + e.z)); });
    return poly;
}
vector<TriangleAttribute> split_polygon(const Polygon& poly) {
    vector<TriangleAttribute> res;
    int n = poly.vertices.size() - 1;
    for(int k = 0; k <= n - 2; k++) {
        TriangleAttribute tri;
        if(k % 2 == 0) {
            tri.vertices[0] = poly.vertices[k / 2];
            tri.vertices[1] = poly.vertices[n - k / 2];
            tri.vertices[2] = poly.vertices[(k + 2) / 2];
        } else {
            tri.vertices[0] = poly.vertices[n - (k - 1) / 2];
            tri.vertices[1] = poly.vertices[(k + 1) / 2];
            tri.vertices[2] = poly.vertices[n - (k + 1) / 2];
        }
        res.push_back(tri);
    }
    return res;
}
float line_dist2D(const vec2& p, const vec2& a, const vec2& b) {
    vec2 n = vec2(a.y - b.y, b.x - a.x);    // xAB
    float d = a.x * b.y - a.y * b.x;
    return (glm::dot(n, p) + d) / glm::length(n);
}
float cal_bary2D(const vec2& a, const vec2& b, const vec2& c, const vec2& p) { return line_dist2D(p, b, c) / line_dist2D(a, b, c); }
vec2 barycentric2D(const vec2& weights, const vec2 vecs[]) { return weights[0] * vecs[0] + weights[1] + vecs[1]; }
vec3 cal_bary3D(const vec3& v, const vec3 vecs[]) {
    vec3 e1 = vecs[2] - vecs[1];
    vec3 e2 = vecs[0] - vecs[2];
    vec3 e3 = vecs[1] - vecs[0];
    vec3 n = glm::cross(e1, e2);
    vec3 d1 = v - vecs[0];
    vec3 d2 = v - vecs[1];
    vec3 d3 = v - vecs[2];
    vec3 weights;
    float area = glm::dot(n, n);
    weights[0] = glm::dot(glm::cross(e1, d3), n) / area;
    weights[1] = glm::dot(glm::cross(e2, d1), n) / area;
    weights[2] = glm::dot(glm::cross(e3, d2), n) / area;
    return weights;
}
vec3 barycentric3D(const vec3& weights, const vec3 vecs[]) { return weights[0] * vecs[0] + weights[1] * vecs[1] + weights[2] * vecs[2]; }
float barycentric3D(const vec3& weights, const float fs[]) {return weights[0] * fs[0] + weights[1] * fs[1] + weights[2] * fs[2]; }
int gamma_encode(float radiance) { return (int)(pow(std::min(1.0f, std::max(0.0f, radiance)), 1.0f/2.2f) * 255.0f); }
vec3 shade(const Scene& scene, const BSDF& bsdf, 
           const vec3& pi, const vec3& normal, const vec3& w_o) {
    vec3 radiance = vec3(0.0f, 0.0f, 0.0f);
    for(int i = 0; i < scene.lights.size(); i++) {
        const Light& light = scene.lights[i];
        vec3 w_i = light.org - pi;
        float dist2light = glm::length(w_i);
        w_i = w_i / dist2light;
        // diffuse: k_diff * c_diff * c_light * 1/r^2 * max(0, dot(n, w_i))
        float k_diff = 5.0f;
        radiance += k_diff * bsdf.diff_color * light.color * 1.0f / (dist2light*dist2light)
                    * std::max(0.0f, glm::dot(normal, w_i));
        // specular: k_spec * c_spec * c_light * 1/r^2 * max(0,dot(h, n))^spec_gloss
        float k_spec = 2.0f;
        vec3 h = glm::normalize(w_o + w_i);
        radiance += k_spec * bsdf.spec_color * light.color * 1.0f / (dist2light*dist2light)
                    * pow(std::max(0.0f, glm::dot(h, normal)), bsdf.spec_gloss);
    }
    return radiance;
}
void screen_mapping(const TriangleAttribute& tri, int w, int h, glm::i32vec2 scr_coords[]) {
    for(int i = 0; i < 3; i++) {
        const vec4& ndc_coord = tri.vertices[i];
        float y = (ndc_coord.y - (-1)) / 2.0f * (h - 1);
        float x = (ndc_coord.x - (-1)) / 2.0f * (w - 1);     
        scr_coords[i] = glm::i32vec2(x, y);
    }
}
void compute_bounding_box(const glm::i32vec2 scr_coords[], int w, int h, int& x_min, int& y_min, int& x_max, int& y_max) {
    x_min = w - 1;  x_max = 0;
    y_min = h - 1;  y_max = 0; 
    for(int i = 0; i < 3; i++) {
        int x = scr_coords[i].x, y = scr_coords[i].y;
        x_min = glm::min(x, x_min);
        x_max = glm::max(x, x_max);
        y_min = glm::min(y, y_min);
        y_max = glm::max(y, y_max); 
    }
    x_min = glm::max(    0, x_min);
    x_max = glm::min(w - 1, x_max);
    y_min = glm::max(    0, y_min);
    y_max = glm::min(h - 1, y_max);     
}
vec3 homo_to_3D(const vec4& hv, const Camera& camera) {
    float tf2 = glm::tan(glm::radians(camera.deg_fov / 2));
    return vec3(hv.x * tf2 * camera.aspect, hv.y * tf2, -hv.w);
}
vector<TriangleAttribute> geometry_process(const Scene& scene, const Camera& camera) {
    vector<TriangleAttribute> tri_attrs;
    mat4 projection = glm::perspective(glm::radians(camera.deg_fov), camera.aspect, camera.z_near, camera.z_far);
    for(int i = 0; i < scene.triangles.size(); i++) {
        const Triangle& tri = scene.triangles[i];
        Polygon poly(tri);
        poly.apply(projection); // geometry transform
        poly = homo_clip(poly); // homogeneous clip
        vector<TriangleAttribute> splited_tris = split_polygon(poly);    // generate new triangles
        for(int j = 0; j < splited_tris.size(); j++) {
            TriangleAttribute& tri_attr = splited_tris[j];
            for(int k = 0; k < 3; k++) {    // derive attributes from old triangles
                tri_attr.wd_pos [k] = homo_to_3D(tri_attr.vertices[k], camera);
                tri_attr.normals[k] = barycentric3D(cal_bary3D(tri_attr.wd_pos[k], tri.vertices), tri.normals);
                tri_attr.bsdf       = tri.bsdf;
            }
        }
        tri_attrs.insert(tri_attrs.end(), splited_tris.begin(), splited_tris.end());
    }
    return tri_attrs;
}
void perspective_divide(TriangleAttribute& tri_attr) { tri_attr.vertices[0] /= tri_attr.vertices[0].w; tri_attr.vertices[1] /= tri_attr.vertices[1].w; tri_attr.vertices[2] /= tri_attr.vertices[2].w; }
void compute_fragments(Image& image, vector<TriangleAttribute>& tri_attrs, Fragment frags[]) {
    int w = image.w, h = image.h;
    for(int i = 0; i < tri_attrs.size(); i++) {
        TriangleAttribute& tri = tri_attrs[i];
        perspective_divide(tri);               // perspective divide
        
        glm::i32vec2 scr_coords[3];
        screen_mapping(tri, w, h, scr_coords);  // screen mapping

        int x_min, y_min, x_max, y_max;
        compute_bounding_box(scr_coords, w, h, x_min, y_min, x_max, y_max); // compute bounding box

        float vertexW [3];
        vec3  vertexPw[3], vertexNw[3];
        for(int j = 0; j < 3; j++) {        // pre-compute for perspective-correct intersection
            vertexW [j] = -1.0f / tri.wd_pos[j].z;
            vertexPw[j] = tri.wd_pos[j] * vertexW[j];
            vertexNw[j] = tri.normals[j] * vertexW[j];
        }
        for(int y = y_min; y <= y_max; y++) {
            for(int x = x_min; x <= x_max; x++) {
                vec3 weights = vec3(
                    cal_bary2D(scr_coords[0], scr_coords[1], scr_coords[2], vec2(x, y)),
                    cal_bary2D(scr_coords[1], scr_coords[2], scr_coords[0], vec2(x, y)),
                    cal_bary2D(scr_coords[2], scr_coords[0], scr_coords[1], vec2(x, y))
                );
                if(weights[0] >= 0 && weights[0] <= 1 &&
                   weights[1] >= 0 && weights[1] <= 1 &&
                   weights[2] >= 0 && weights[2] <= 1   ) {
                    float vertexZ[3] = {tri.vertices[0].z, tri.vertices[1].z, tri.vertices[2].z};
                    float perspective_z = barycentric3D(weights, vertexZ);
                    if(perspective_z < image.z_buffer[image.get_index(x, y)]) {
                        image.z_buffer[image.get_index(x, y)] = perspective_z;
                        Fragment frag;
                        frag.wd_pos = barycentric3D(weights, vertexPw) / barycentric3D(weights, vertexW);
                        frag.normal = glm::normalize(barycentric3D(weights, vertexNw) / barycentric3D(weights, vertexW));
                        frag.bsdf   = tri.bsdf;
                        frags[image.get_index(x, y)] = frag;
                    }
                }                
            }
        }
    }
}
void rasterization(Image& image, const Scene& scene, const Camera& camera) {
    int w = image.w, h = image.h;
    vector<TriangleAttribute> tri_attrs = geometry_process(scene, camera);  // after homo clip, before perspective div
    Fragment* frags = new Fragment[w*h];
    compute_fragments(image, tri_attrs, frags);
    {
        fixed_thread_pool pool(32);
        for(int i = 0; i < w*h; i++) {
            pool.execute([i, &image, &frags, &scene](){
                if(image.z_buffer[i] != std::numeric_limits<float>::infinity()) {
                    const Fragment& frag = frags[i];
                    int x, y;
                    vec3 radiance = shade(scene, frag.bsdf, frag.wd_pos, frag.normal, -glm::normalize(frag.wd_pos));
                    u8vec3 col(gamma_encode(radiance.r), gamma_encode(radiance.g), gamma_encode(radiance.b));
                    image.get_xy(i, x, y);
                    image.fill_color(x, y, col);
                    
                }
            });
        }
    }
    delete[] frags;
    image.write("rendering_pipeline.ppm");
}
void init_scene(Scene& scene) {
    BSDF bsdf;
    bsdf.diff_color = vec3(0.05f, 0.4f, 0.1f);
    bsdf.spec_color = vec3( 1.0f, 1.0f, 1.0f);
    bsdf.spec_gloss = 320.0f;

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
    Camera camera(90.0f, 0.3f, 100.0f, (float)w / h);
    Scene scene;
    init_scene(scene);
    rasterization(image, scene, camera);
    return 0;
}
