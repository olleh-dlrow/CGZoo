/*
c++ homo_clip.cpp -std=c++11 -o homo_clip && ./homo_clip
*/
#include<glm/glm.hpp>
#include<glm/ext.hpp>
#include<iostream>
#include<vector>
#include<functional>
using std::function;
using std::vector;
using glm::vec3;
using glm::vec4;
using glm::mat4;
void print(const glm::vec4& vec) {
    for(int i = 0; i < 4; i++) {
        printf("%7.3f,", vec[i]);
    }std::cout << "\n";
}
class Camera {
public:
    Camera(float fov, float zn, float zf, float a):deg_fov(fov), z_near(zn), z_far(zf), aspect(a) {}
    float deg_fov;  // degree
    float z_near;
    float z_far;
    float aspect;
};
class Polygon {
public:
    vector<vec4> vertices;
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
        
        if(!in_same_side(s_sign, inner_sign) && !in_same_side(e_sign, inner_sign)) {
            // 1. s out, e out: output None
            continue;
        } else if(in_same_side(s_sign, inner_sign) && in_same_side(e_sign, inner_sign)) {
            // 2. s in, e in: output e
            res.vertices.push_back(e);
        } else if(!in_same_side(s_sign, inner_sign) && in_same_side(e_sign, inner_sign)) {
            // 3. s out, e in: output i, e
            float t = cal_t(s, e);
            vec4 intersect = (1 - t) * s + t * e;
            res.vertices.push_back(intersect);
            res.vertices.push_back(e);
        } else {
            // 4. s in, e out: output i
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

int main() {
    Camera camera(90.0f, 1.0f, 100.0f, 1.0f);
    mat4 projection = glm::perspective(glm::radians(camera.deg_fov), camera.aspect, camera.z_near, camera.z_far);
    vec4 v1 = vec4( 0, 0, -4, 1);
    vec4 v2 = vec4( 5, 0,  1, 1);
    vec4 v3 = vec4(-5, 0,  1, 1);
    Polygon poly;
    poly.vertices.push_back(projection * v1);
    poly.vertices.push_back(projection * v2);
    poly.vertices.push_back(projection * v3);
    printf("homogeneous coordinates before clip:\n");
    print(poly.vertices[0]);
    print(poly.vertices[1]);
    print(poly.vertices[2]);
    printf("\n");

    poly = homo_clip(poly);    
    printf("NDC coordinates after clip:\n");
    for(int i = 0; i < poly.vertices.size(); i++) {
        const vec4& v = poly.vertices[i];
        print(v/v.w);
    }
    return 0;
}
