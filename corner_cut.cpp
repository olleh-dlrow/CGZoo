/*
clang++ corner_cut.cpp -o corner_cut && ./corner_cut
*/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<vector>
using std::vector;

class Vec2;
Vec2 lerp(const Vec2& v1, const Vec2& v2, float t);

class Color{
public:
    Color():r(255), g(255), b(255) {}
    Color(unsigned char _r, unsigned char _g, unsigned char _b):r(_r), g(_g), b(_b) {}
    void draw(unsigned char _r, unsigned char _g, unsigned char _b) {
        r = _r; g = _g; b = _b;
    }
    unsigned char r;
    unsigned char g;
    unsigned char b;
};

class Vec2{
public:
    Vec2():x(0), y(0){}
    Vec2(float _x, float _y):x(_x), y(_y) {}
    Vec2 operator+(const Vec2& v) const {return Vec2(v.x + x, v.y + y);}
    Vec2 operator-(const Vec2& v) const {return Vec2(x - v.x, y - v.y);}
    Vec2 operator*(float t) const {return Vec2(x * t, y * t);}
    float x;
    float y;
};

class Line{
public:
    Line(const Vec2& _v0, const Vec2& _v1) {
        v0 = _v0.x <= _v1.x ? _v0 : _v1;
        v1 = _v0.x > _v1.x ? _v0 : _v1;
    }
    inline float f(const Vec2& v) const {
        // (y0 - y1)*x + (x1 - x0)*y + x0y1 - x1y0 = 0
        return (v0.y - v1.y) * v.x + (v1.x - v0.x) * v.y + v0.x * v1.y - v1.x * v0.y;
    }
    void print() const {
        printf("(%f,%f)<->(%f,%f)\n", v0.x, v0.y, v1.x, v1.y);
    }
    Vec2 v0;
    Vec2 v1;
};

class Polygon {
public:
    Polygon() {}
    Polygon corner_cut(int nb_cut) {
        assert(nb_cut > 0);
        Polygon poly;
        int sz = vertices.size();
        for(int s = 0; s < sz; s++) {
            int e = (s + 1) % sz;
            float t = 1.0 / (nb_cut + 1);
            poly.vertices.push_back(lerp(vertices[s], vertices[e], t));
            if(nb_cut != 1)
                poly.vertices.push_back(lerp(vertices[e], vertices[s], t));
        }
        return poly;
    }
    vector<Vec2> vertices;
};

class Image {
public:
    Image(int _w, int _h):w(_w), h(_h) {
        colors = new Color[w*h];
    }

    void fill_color(int x, int y, const Color& color) {
        int index = (h - 1 - y) * w + x;
        if(x >= 0 && y >= 0 && x < w && y < h)
            colors[index] = color;
    }

    void draw_polygon(const Polygon& polygon, const Color& color) {
        int nb_vertices = polygon.vertices.size();
        for(int s = 0; s < nb_vertices; s++) {
            int e = (s + 1) % nb_vertices;
            Line line(polygon.vertices[s], polygon.vertices[e]);
            draw_line(line, color);
        }
    }

    void draw_line(const Line& line, const Color& color) {
        Vec2 v0 = line.v0;
        Vec2 v1 = line.v1;

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
                    Vec2 vec(x + 0.5 * sign, y + 1);
                    if(line.f(vec) > 0) {
                        x += sign;
                    }
                }
            } else {    // -1 <= k <= 1
                int y = v0.y;
                for (int x = v0.x; x <= v1.x; x++)
                {
                    fill_color(x, y, color);
                    Vec2 vec(x + 1, y + 0.5 * sign);
                    if(line.f(vec) * sign < 0) {
                        y += sign;
                    }
                }
            }
        }
    }
    void write(const char* filename) {
        FILE* f = fopen(filename, "w");
        fprintf(f, "P3\n%d %d\n255\n", w, h);
        for (int i = 0; i < w*h; i++) {
            fprintf(f, "%d %d %d ", colors[i].r, colors[i].g, colors[i].b);
        }
        fclose(f);  
    }
    ~Image() {
        delete[] colors;
    }
    Color* colors;
    int w;
    int h;
};

Vec2 lerp(const Vec2& v1, const Vec2& v2, float t) {
    // v = v1 + (v2 - v1) * t
    return v1 + (v2 - v1) * t;
}

float deg2rad(float deg) { return deg * M_PI / 180.0; }

Polygon get_triangle(const Vec2& center, float r) {
    Polygon polygon;
    polygon.vertices.push_back(Vec2(center.x + r * cos(deg2rad(90)), center.y + r * sin(deg2rad(90))));
    polygon.vertices.push_back(Vec2(center.x + r * cos(deg2rad(210)), center.y + r * sin(deg2rad(210))));
    polygon.vertices.push_back(Vec2(center.x + r * cos(deg2rad(-30)), center.y + r * sin(deg2rad(-30))));
    return polygon;
}

Polygon get_square(const Vec2& center, float a) {
    Polygon polygon;
    polygon.vertices.push_back(Vec2(center.x - a/2, center.y - a/2));
    polygon.vertices.push_back(Vec2(center.x + a/2, center.y - a/2));
    polygon.vertices.push_back(Vec2(center.x + a/2, center.y + a/2));
    polygon.vertices.push_back(Vec2(center.x - a/2, center.y + a/2));
    return polygon;
}

int main(int argc, char* argv[])
{
    int w = 1024, h = 1024;
    Vec2 center(w / 2, h / 2);
    Image image(w, h);
    Color line_color(0, 0, 0);
    float r = 300;

    Polygon polygon = get_square(center, 600);
    image.draw_polygon(polygon, line_color);
    int cnt_iter = 10;
    for(int i = 0; i < cnt_iter; i++) {
        polygon = polygon.corner_cut(1);
        image.draw_polygon(polygon, line_color);
    }
    image.write("corner_cut.ppm");
    return 0;
}