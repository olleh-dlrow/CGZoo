#include<iostream>
#include<vector>
#include<string>
#include<map>
#include<set>
#include<memory>
#include<unordered_map>
#include<array>
#include<fstream>
#include<random>

using namespace std;

using Float = float;
struct Vec2 {
    Float x;
    Float y;

    Float operator[](int index) {
        switch(index) {
            case 0:
                return x;
            case 1:
                return y;
            default:
                assert(0);
        }
        return 0;
    }

    Vec2 operator+(const Vec2& rhs) const {
        return {x + rhs.x, y + rhs.y};
    }
};

Vec2 operator*(Float lhs, const Vec2& rhs) {
    return {lhs * rhs.x, lhs * rhs.y};
}

array<float, 3> SampleUniformTriangle(Vec2 u) {
    float b0, b1;
    if (u[0] < u[1]) {
        b0 = u[0] / 2;
        b1 = u[1] - b0;
    } else {
        b1 = u[1] / 2;
        b0 = u[0] - b1;
    }
    return {b0, b1, 1 - b0 - b1};
}

// ERROR!!! not uniform
array<float, 3> SampleUniformTriangle2(Vec2 u) {
    return {1 - u.x, u.x * (1 - u.y), u.x * u.y};
}

array<float, 3> SampleUniformTriangle3(Vec2 u) {
    if(u.x + u.y > 1) {
        u.x = 1 - u.x;
        u.y = 1 - u.y;
    }
    return {u.x, u.y, 1 - u.x - u.y};
}

struct Vec3 {
    Float x = 0;
    Float y = 0;
    Float z = 0;

    Vec3() {

    }

    Vec3(Float _x, Float _y, Float _z):x(_x), y(_y), z(_z) {

    }

    Float operator[](int index) const {
        switch(index) {
            case 0:
                return x;
            case 1:
                return y;
            case 2:
                return z;
        }
        return -1;
    }

    Vec3 cross(const Vec3& rhs) const {
        return Vec3{
            y * rhs.z - z * rhs.y,
            z * rhs.x - x * rhs.z,
            x * rhs.y - y * rhs.x
        };
    }

    Float dot(const Vec3& rhs) const {
        return x * rhs.x + y * rhs.y + z * rhs.z;
    }

    Vec3 operator*(Float rhs) const {
        return Vec3{x * rhs, y * rhs, z * rhs};
    }

    Vec3 operator+(const Vec3& rhs) const {
        return Vec3{x + rhs.x, y + rhs.y, z + rhs.z};
    }

    Vec3 operator-(const Vec3& rhs) const {
        return Vec3{x - rhs.x, y - rhs.y, z - rhs.z};
    }

    Vec3& operator+=(const Vec3& rhs) {
        *this = *this + rhs;
        return *this;
    }

    Vec3 operator/(Float rhs) const {
        return Vec3{x/rhs, y/rhs, z/rhs};
    }

    std::string to_string() const {
        return "(" + std::to_string(x) + ", " + std::to_string(y) + ", "
            + std::to_string(z) + ")";
    }

    static Vec3 min(const Vec3& lhs, const Vec3& rhs) {
        return {
            std::min(lhs.x, rhs.x),
            std::min(lhs.y, rhs.y),
            std::min(lhs.z, rhs.z),
        };
    }

    static Vec3 max(const Vec3& lhs, const Vec3& rhs) {
            return {
                std::max(lhs.x, rhs.x),
                std::max(lhs.y, rhs.y),
                std::max(lhs.z, rhs.z),
            };
    }
};

Vec3 operator*(Float v, const Vec3& rhs) {
    return rhs * v;
}

Vec3 operator/(Float v, const Vec3& rhs) {
    return Vec3(v/rhs.x, v/rhs.y, v/rhs.z);
}


struct Image {
    Image(int w, int h):
        width(w), height(h) {
        colors = std::vector<Vec3>(w * h, {0.0f, 0.0f, 0.0f});
    }

    void set_pixel(int x, int y, Vec3 col) {
        int idx = (height - 1 - y) * width + x;
        colors[idx] = col;
    }

    void write_ppm(const std::string& filename) {
        std::ofstream file(filename);
        file << "P3\n" << width << " " << height << "\n255\n";
        for(int i = 0; i < width * height; i++) {
            file << colors[i].x << " " << colors[i].y << " " 
                << colors[i].z << " ";
        }
    }
    std::vector<Vec3> colors;
    int width, height;
};

Float get_uniform01()  {
    random_device dev;
    std::mt19937 rng(dev());
    std::uniform_real_distribution<float> dist(0.f, 1.f);
    return dist(rng);
}

int main() {
    Image img(1024, 1024);
    Vec2 v[3] = {
        {100, 200},
        {400, 900},
        {800, 540}
    };
    int cnt = 100000;
    while(cnt--) {
        auto b = SampleUniformTriangle3({get_uniform01(), get_uniform01()});
        Vec2 p = b[0] * v[0] + b[1] * v[1] + b[2] * v[2];
        img.set_pixel(p.x, p.y, {255, 255, 255});
    }
    img.write_ppm("triangle_sample.ppm");

    return 0;
}
