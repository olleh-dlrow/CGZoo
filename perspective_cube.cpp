/*
clang++ perspective_cube.cpp -o perspective_cube && ./perspective_cube
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
    Vec2(float _x, float _y):x(_x), y(_y) {}
    float x;
    float y;
};

class Vec3: public Vec2{
public:
    Vec3():Vec2(), z(0) {}
    Vec3(float _x, float _y, float _z):Vec2(_x, _y), z(_z) {}
    float z;
};

class Cube {
public:
    Cube() {}

    Vec3 vertices[8];
    int edges[12][2];
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

// map [0,1]x[0,1] to [0,w]x[0,h]
Vec2 map(const Vec2& v, int w, int h) {
    return Vec2(round(v.x * w - 0.5), round(v.y * h - 0.5));
}

void init_cube(Cube& cube, int x_min, int x_max, int y_min, int y_max, int z_min, int z_max) {
    cube.vertices[0].x=x_min;cube.vertices[0].y=y_min;cube.vertices[0].z=z_min;
    cube.vertices[1].x=x_max;cube.vertices[1].y=y_min;cube.vertices[1].z=z_min;
    cube.vertices[2].x=x_min;cube.vertices[2].y=y_max;cube.vertices[2].z=z_min;
    cube.vertices[3].x=x_max;cube.vertices[3].y=y_max;cube.vertices[3].z=z_min;
    cube.vertices[4].x=x_min;cube.vertices[4].y=y_min;cube.vertices[4].z=z_max;
    cube.vertices[5].x=x_max;cube.vertices[5].y=y_min;cube.vertices[5].z=z_max;
    cube.vertices[6].x=x_min;cube.vertices[6].y=y_max;cube.vertices[6].z=z_max;
    cube.vertices[7].x=x_max;cube.vertices[7].y=y_max;cube.vertices[7].z=z_max;
    cube.edges[0 ][0]=0;cube.edges[0 ][1]=1;
    cube.edges[1 ][0]=1;cube.edges[1 ][1]=3;
    cube.edges[2 ][0]=2;cube.edges[2 ][1]=3;
    cube.edges[3 ][0]=0;cube.edges[3 ][1]=2;
    cube.edges[4 ][0]=4;cube.edges[4 ][1]=5;
    cube.edges[5 ][0]=5;cube.edges[5 ][1]=7;
    cube.edges[6 ][0]=6;cube.edges[6 ][1]=7;
    cube.edges[7 ][0]=4;cube.edges[7 ][1]=6;
    cube.edges[8 ][0]=0;cube.edges[8 ][1]=4;
    cube.edges[9 ][0]=1;cube.edges[9 ][1]=5;
    cube.edges[10][0]=3;cube.edges[10][1]=7;
    cube.edges[11][0]=2;cube.edges[11][1]=6;
}

int main(int argc, char* argv[])
{
    int w = 1024, h = 1024;
    float x_min = -300, x_max = 300;
    float y_min = -300, y_max = 300;
    float z_min =  500, z_max = 1100;

    float frame_width = 600, frame_height = 600;
    float frame_left   = -frame_width/2,  frame_right = frame_width/2;
    float frame_bottom = -frame_height/2, frame_top   = frame_height/2;
    float frame_near = 120;
    float frame_far = 700;

    Vec3 eye(0, 0, 0);
    Image image(w, h);
    Color line_color(0, 0, 0);

    float width = 800, height = 800, depth = 800;
    float offset = 1500, offz = 2000;
    Vec3 centers[] = {
        Vec3(0, 0,              offz),
        Vec3(offset, 0,         offz),
        Vec3(-offset, 0,        offz),
        Vec3(0, offset,         offz),
        Vec3(0, -offset,        offz),
        Vec3(offset, offset,    offz),
        Vec3(-offset, -offset,  offz),
        Vec3(-offset, offset,   offz),
        Vec3(offset, -offset,   offz)
    };
    int len = sizeof(centers) / sizeof(Vec3);

    for(int l = 0; l < len; l++) {
        x_min = centers[l].x - width/2;  x_max = centers[l].x + width/2;
        y_min = centers[l].y - height/2; y_max = centers[l].y + height/2;
        z_min = centers[l].z - depth/2;  z_max = centers[l].z + depth/2;
        Cube cube;
        init_cube(cube, x_min, x_max, y_min, y_max, z_min, z_max);

        Vec3 per_vertices[8];
        for(int i = 0; i < 8; i++) {
            per_vertices[i].x = (cube.vertices[i].x) * (frame_near) / (cube.vertices[i].z);
            per_vertices[i].x = (per_vertices[i].x - frame_left) / (frame_right - frame_left);

            per_vertices[i].y = (cube.vertices[i].y) * (frame_near) / (cube.vertices[i].z);
            per_vertices[i].y = (per_vertices[i].y - frame_bottom) / (frame_top - frame_bottom);

            per_vertices[i].z = frame_near + frame_far - (frame_near * frame_far) / per_vertices[i].z;
        }

        for(int i = 0; i < 12; i++) {
            Vec2 s = per_vertices[cube.edges[i][0]];
            Vec2 e = per_vertices[cube.edges[i][1]];
            Vec2 s_px = map(s, w, h);
            Vec2 e_px = map(e, w, h);
            Line line(s_px, e_px);
            image.draw_line(line, line_color);
        }
         
    }
   
    image.write("perspective_cube.ppm");
    return 0;
}