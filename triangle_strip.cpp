/*
c++ triangle_strip.cpp -std=c++11 -o triangle_strip && ./triangle_strip
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
vector<Polygon> split_polygon(const Polygon& poly) {
    vector<Polygon> res;
    int n = poly.vertices.size() - 1;
    for(int k = 0; k <= n - 2; k++) {
        Polygon tri;
        if(k % 2 == 0) {
            tri.vertices.push_back(poly.vertices[k / 2]);
            tri.vertices.push_back(poly.vertices[n - k / 2]);
            tri.vertices.push_back(poly.vertices[(k + 2) / 2]);
        } else {
            tri.vertices.push_back(poly.vertices[n - (k - 1) / 2]);
            tri.vertices.push_back(poly.vertices[(k + 1) / 2]);
            tri.vertices.push_back(poly.vertices[n - (k + 1) / 2]);
        }
        res.push_back(tri);
    }
    return res;
}
int main() {
    int w = 1024, h = 1024;
    Image image(w, h);
    u8vec3 black(0, 0, 0);
    Polygon poly;
    poly.vertices.push_back(vec2(512, 900));
    poly.vertices.push_back(vec2(600, 910));
    poly.vertices.push_back(vec2(950, 512));
    poly.vertices.push_back(vec2(912, 150));
    poly.vertices.push_back(vec2(200, 120));
    poly.vertices.push_back(vec2(150, 150));
    poly.vertices.push_back(vec2( 50, 400));
    // image.draw_polygon(poly, black);
    vector<Polygon> triangles = split_polygon(poly);
    for(int i = 0; i < triangles.size(); i++)image.draw_polygon(triangles[i], black);
    image.write("triangle_strip.ppm");
    return 0;
}