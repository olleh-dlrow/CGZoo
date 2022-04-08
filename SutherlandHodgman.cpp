/*
c++ SutherlandHodgman.cpp -std=c++11 -o SutherlandHodgman && ./SutherlandHodgman
*/
#include<stdio.h>
#include<vector>
#include<functional>
#include<glm/ext.hpp>
#include<glm/glm.hpp>
using std::vector;
using std::function;
using glm::vec2;
using glm::u8vec3;
class Line{
public:
    Line(const vec2& _v0, const vec2& _v1) {
        v0 = _v0.x <= _v1.x ? _v0 : _v1;
        v1 = _v0.x > _v1.x ? _v0 : _v1;
    }
    inline float f(const vec2& v) const {
        // (y0 - y1)*x + (x1 - x0)*y + x0y1 - x1y0 = 0
        return (v0.y - v1.y) * v.x + (v1.x - v0.x) * v.y + v0.x * v1.y - v1.x * v0.y;
    }
    void print() const {
        printf("(%f,%f)<->(%f,%f)\n", v0.x, v0.y, v1.x, v1.y);
    }
    vec2 v0;
    vec2 v1;
};
class Polygon {
public:
    Polygon() {}
    vector<vec2> vertices;
};
class Image {
public:
    Image(int _w, int _h):w(_w), h(_h) { 
        colors = new glm::u8vec3[w*h];
        for(int i = 0; i < w*h; i++) {
            colors[i] = u8vec3(255, 255, 255);
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

    void draw_line(const Line& line, const glm::u8vec3& color) {
        vec2 v0 = line.v0;
        vec2 v1 = line.v1;

        int sign = v1.y - v0.y >= 0 ? 1 : -1;
        // inf
        if(v0.x == v1.x) {
            int lower = v0.y <= v1.y ? v0.y : v1.y;
            int upper = v0.y > v1.y ? v0.y : v1.y;
            int x = v0.x;
            for(int y = lower; y <= upper; y++) {
                fill_color(x, y, color);
            }            
        }
        // k = (y1 - y0) / (x1 - x0) 
        else {
            float k = (v1.y - v0.y) / (v1.x - v0.x);
            // k > 1 or k < 1
            if(k > 1 || k < -1) {
                int x = v0.y < v1.y ? v0.x : v1.x;
                int lower = v0.y < v1.y ? v0.y : v1.y;
                int upper = v0.y > v1.y ? v0.y : v1.y;
                for(int y = lower; y <= upper; y++) {
                    fill_color(x, y, color);
                    vec2 vec(x + 0.5 * sign, y + 1);
                    if(line.f(vec) > 0) {
                        x += sign;
                    }
                }
            } else {    // -1 <= k <= 1
                int y = v0.y;
                for (int x = v0.x; x <= v1.x; x++)
                {
                    fill_color(x, y, color);
                    vec2 vec(x + 1, y + 0.5 * sign);
                    if(line.f(vec) * sign < 0) {
                        y += sign;
                    }
                }
            }
        }
    }
    void draw_polygon(const Polygon& polygon, const u8vec3& color) {
        int nb_vertices = polygon.vertices.size();
        for(int s = 0; s < nb_vertices; s++) {
            int e = (s + 1) % nb_vertices;
            Line line(polygon.vertices[s], polygon.vertices[e]);
            draw_line(line, color);
        }
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
inline bool in_same_side(float s1, float s2){ return s1 * s2 >= 0; }
// input: 待裁剪的多边形
// sign:  判断点的相对位置
// inner: 在裁剪直线内部的点，如果sign(v) == sign(inner)，则该点在内部
// cal_t: 根据s和e计算插值t, 从而得到交点
Polygon sutherland_hodgman(const Polygon& input, 
                           function<float(const vec2& v)> sign, 
                           const vec2& inner,
                           function<float(const vec2& s, const vec2& e)> cal_t) {
    Polygon res;
    for(int i = 0; i < input.vertices.size(); i++) {
        const vec2& s = input.vertices[i];
        const vec2& e = input.vertices[(i + 1) % input.vertices.size()];
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
            vec2 intersect = (1 - t) * s + t * e;
            res.vertices.push_back(intersect);
            res.vertices.push_back(e);
        } else {
            // 4. s in, e out: output i
            float t = cal_t(s, e);
            vec2 intersect = (1 - t) * s + t * e;
            res.vertices.push_back(intersect);
        }
    }
    return res;
}
int main() {
    int w = 1024, h = 1024;
    Image   image(w, h);
    u8vec3  black = u8vec3(0  , 0, 0);
    u8vec3  red   = u8vec3(255, 0, 0);
    Polygon win, poly;
    vec2    inner = vec2(w / 2, h / 2);
    float wx0 = 112, wx1 = 912, wy0 = 212, wy1 = 812;
    win.vertices.push_back(vec2(wx0, wy0));
    win.vertices.push_back(vec2(wx0, wy1));
    win.vertices.push_back(vec2(wx1, wy1));
    win.vertices.push_back(vec2(wx1, wy0));

    poly.vertices.push_back(vec2(512, 900));
    poly.vertices.push_back(vec2(600, 700));
    poly.vertices.push_back(vec2(950, 512));
    poly.vertices.push_back(vec2(912, 150));
    poly.vertices.push_back(vec2(200, 300));
    poly.vertices.push_back(vec2(150, 250));
    poly.vertices.push_back(vec2( 50, 400));
    image.draw_polygon(win , black);
    image.draw_polygon(poly, black);
    
    poly = sutherland_hodgman(poly, 
        [=](const vec2& v)->float{return v.y - wy1;}, inner, 
        [=](const vec2& s, const vec2& e)->float{ return (wy1 - s.y) / (e.y - s.y);});
    poly = sutherland_hodgman(poly, 
        [=](const vec2& v)->float{return v.x - wx1;}, inner, 
        [=](const vec2& s, const vec2& e)->float{ return (wx1 - s.x) / (e.x - s.x);});
    poly = sutherland_hodgman(poly, 
        [=](const vec2& v)->float{return v.y - wy0;}, inner, 
        [=](const vec2& s, const vec2& e)->float{ return (wy0 - s.y) / (e.y - s.y);});
    poly = sutherland_hodgman(poly, 
        [=](const vec2& v)->float{return v.x - wx0;}, inner, 
        [=](const vec2& s, const vec2& e)->float{ return (wx0 - s.x) / (e.x - s.x);});   
    image.draw_polygon(poly, red);
    image.write("SutherlandHodgman.ppm");
    return 0;
}