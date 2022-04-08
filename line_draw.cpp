/*
clang++ line_draw.cpp -o line_draw && ./line_draw
*/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>

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
    Vec2(int _x, int _y):x(_x), y(_y) {}
    int x;
    int y;
};

class Line{
public:
    Line(const Vec2& _v0, const Vec2& _v1) {
        v0 = _v0.x <= _v1.x ? _v0 : _v1;
        v1 = _v0.x > _v1.x ? _v0 : _v1;
    }
    int f(Vec2& v){
        // (y0 - y1)*x + (x1 - x0)*y + x0y1 - x1y0 = 0
        return (v0.y - v1.y) * v.x + (v1.x - v0.x) * v.y + v0.x * v1.y - v1.x * v0.y;
    }
    void print() const {
        printf("(%d,%d)<->(%d,%d)\n", v0.x, v0.y, v1.x, v1.y);
    }
    Vec2 v0;
    Vec2 v1;
};

void draw_line(Color colors[], Line& line,  const Color& color, int h, int w) {
    Vec2 v0 = line.v0, v1 = line.v1;
    int sign = v1.y - v0.y >= 0 ? 1 : -1;
    // inf
    if(v0.x == v1.x) {
        int lower = v0.y <= v1.y ? v0.y : v1.y;
        int upper = v0.y > v1.y ? v0.y : v1.y;
        int x = v0.x;
        for(int y = lower; y <= upper; y++) {
            colors[(h - 1 - y) * w + x] = color;
        }            
    }
    // k = (y1 - y0) / (x1 - x0) 
    else {
        float k = (float)(v1.y - v0.y) / (float)(v1.x - v0.x);
        // k > 1 or k < 1
        if(k > 1 || k < -1) {
            int x = v0.y < v1.y ? v0.x : v1.x;
            int lower = v0.y < v1.y ? v0.y : v1.y;
            int upper = v0.y > v1.y ? v0.y : v1.y;
            for(int y = lower; y <= upper; y++) {
                colors[(h - 1 - y) * w + x] = color;
                Vec2 vec(x + 0.5 * sign, y + 1);
                if(line.f(vec) > 0) {
                    x += sign;
                }
            }
        } else {    // -1 <= k <= 1
            int y = v0.y;
            for (int x = v0.x; x <= v1.x; x++)
            {
                colors[(h - 1 - y) * w + x] = color;
                Vec2 vec(x + 1, y + 0.5 * sign);
                if(line.f(vec) * sign < 0) {
                    y += sign;
                }
            }
        }
    }

}



int main(int argc, char* argv[])
{
    int w = 1024, h = 1024;
    Color* colors = new Color[w*h];
    for(int i = 0; i < w * h; i++) {
        colors[i].draw(0, 0, 0);
    }

    Vec2 center(w / 2, h / 2);
    int r = 300;
    int theta = 0;
    int delta_theta = 5;
    Color select_colors[] = {
        Color(255, 0, 0),
        Color(0, 255, 0),
        Color(0, 0, 255),
        Color(255, 255, 0),
        Color(255, 0, 255),
        Color(0, 255, 255),
        Color(255, 255, 255)
    };
    int iter_color = 0;
    int len = sizeof(select_colors) / sizeof(Color);

    while(theta < 360) {
        Vec2 p(center.x + r * cos(theta), center.y + r * sin(theta));
        float rad = theta * M_PI / 180.0;
        Line line(center, Vec2(center.x + r * cos(rad), center.y + r * sin(rad)));
        draw_line(colors, line, select_colors[iter_color], h, w);
        iter_color = (iter_color + 1) % len;
        theta += delta_theta;
    }

    FILE* f = fopen("line_draw.ppm", "w");
    fprintf(f, "P3\n%d %d\n255\n", w, h);
    for (int i = 0; i < w*h; i++) {
        fprintf(f, "%d %d %d ", colors[i].r, colors[i].g, colors[i].b);
    }
    return 0;
}