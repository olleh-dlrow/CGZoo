/*
c++ sobol_sample.cpp -O3 -std=c++20 -o sobol_sample && ./sobol_sample
*/

#include<iostream>
#include<unordered_map>
#include<vector>
#include<array>
#include<string>
#include<fstream>
#include<sstream>
#include<queue>
#include<memory>
#include<cmath>
#include<random>
#include<numeric>
#include<chrono>
#include<thread>
#include<cassert>
#include<algorithm>

using Float = float;

struct SobolSampler {
    SobolSampler() {

    }

    Float sample(unsigned N, unsigned D) const {
        return P[N][D];
    }

    static double **sobol_points(unsigned N, unsigned D, const char *dir_file)
    {
        std::ifstream infile(dir_file,std::ios::in);
        if (!infile) {
            std::cout << "Input file containing direction numbers cannot be found!\n";
            exit(1);
        }
        char buffer[1000];
        infile.getline(buffer,1000,'\n');
        
        // L = max number of bits needed 
        unsigned L = (unsigned)ceil(log((double)N)/log(2.0)); 

        // C[i] = index from the right of the first zero bit of i
        unsigned *C = new unsigned [N];
        C[0] = 1;
        for (unsigned i=1;i<=N-1;i++) {
            C[i] = 1;
            unsigned value = i;
            while (value & 1) {
            value >>= 1;
            C[i]++;
            }
        }
        
        // POINTS[i][j] = the jth component of the ith point
        //                with i indexed from 0 to N-1 and j indexed from 0 to D-1
        double **POINTS = new double * [N];
        for (unsigned i=0;i<=N-1;i++) POINTS[i] = new double [D];
        for (unsigned j=0;j<=D-1;j++) POINTS[0][j] = 0; 

        // ----- Compute the first dimension -----
        
        // Compute direction numbers V[1] to V[L], scaled by pow(2,32)
        unsigned *V = new unsigned [L+1]; 
        for (unsigned i=1;i<=L;i++) V[i] = 1 << (32-i); // all m's = 1

        // Evalulate X[0] to X[N-1], scaled by pow(2,32)
        unsigned *X = new unsigned [N];
        X[0] = 0;
        for (unsigned i=1;i<=N-1;i++) {
            X[i] = X[i-1] ^ V[C[i-1]];
            POINTS[i][0] = (double)X[i]/pow(2.0,32); // *** the actual points
            //        ^ 0 for first dimension
        }
        
        // Clean up
        delete [] V;
        delete [] X;
        
        
        // ----- Compute the remaining dimensions -----
        for (unsigned j=1;j<=D-1;j++) {
            
            // Read in parameters from file 
            unsigned d, s;
            unsigned a;
            infile >> d >> s >> a;
            unsigned *m = new unsigned [s+1];
            for (unsigned i=1;i<=s;i++) infile >> m[i];
            
            // Compute direction numbers V[1] to V[L], scaled by pow(2,32)
            unsigned *V = new unsigned [L+1];
            if (L <= s) {
            for (unsigned i=1;i<=L;i++) V[i] = m[i] << (32-i); 
            }
            else {
            for (unsigned i=1;i<=s;i++) V[i] = m[i] << (32-i); 
            for (unsigned i=s+1;i<=L;i++) {
            V[i] = V[i-s] ^ (V[i-s] >> s); 
            for (unsigned k=1;k<=s-1;k++) 
            V[i] ^= (((a >> (s-1-k)) & 1) * V[i-k]); 
            }
            }
            
            // Evalulate X[0] to X[N-1], scaled by pow(2,32)
            unsigned *X = new unsigned [N];
            X[0] = 0;
            for (unsigned i=1;i<=N-1;i++) {
            X[i] = X[i-1] ^ V[C[i-1]];
            POINTS[i][j] = (double)X[i]/pow(2.0,32); // *** the actual points
            //        ^ j for dimension (j+1)
        }
            
            // Clean up
            delete [] m;
            delete [] V;
            delete [] X;
        }
        delete [] C;
        
        return POINTS;
    }

    static double** P;
    const static unsigned MAXN = 1024;
    const static unsigned MAXD = 7;
};

double** SobolSampler::P = SobolSampler::sobol_points(
    SobolSampler::MAXN, SobolSampler::MAXD, "new-joe-kuo-6.21201");

struct Vec3 {
    Float x = 0;
    Float y = 0;
    Float z = 0;

    Vec3() {

    }

    Vec3(Float v):
        x(v), y(v), z(v) {

    }

    Vec3(Float _x, Float _y, Float _z):x(_x), y(_y), z(_z) {

    }
    Float magnitude() const {
        return std::sqrt(x*x + y*y + z*z);
    }

    Float sqr_magnitude() const {
        return x*x + y*y + z*z;
    }

    Vec3 operator*(const Vec3& rhs) const {
        return Vec3(
            x * rhs.x,
            y * rhs.y,
            z * rhs.z
        );
    }

    Float operator[](int index) const {
        switch(index) {
            case 0:
                return x;
            case 1:
                return y;
            case 2:
                return z;
            default:
                assert(0);
        }
        return 0;
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

int gamma_encode(Float radiance, bool gamma) {
    if(gamma)
        return (int)(pow(std::min(1.0f, std::max(0.0f, radiance)), 1.0f/2.2f) * 255.0f);
    else 
        return (int)(radiance * 255.0f);
}

struct Image {
    Image(int w, int h):
        width(w), height(h) {
        colors = std::vector<Vec3>(w * h, {0, 0, 0});
    }

    int idx(int x, int y) const {
        return (height - 1 - y) * width + x;
    }

    Vec3& operator[](int index) {
        return colors[index];
    }

    void write_ppm(const std::string& filename, bool gamma) {
        std::ofstream file(filename);
        file << "P3\n" << width << " " << height << "\n255\n";
        for(int i = 0; i < width * height; i++) {
            int col[3] = {
                gamma_encode(colors[i].x, gamma),
                gamma_encode(colors[i].y, gamma),
                gamma_encode(colors[i].z, gamma),
            };
            file << col[0] << " " << col[1] << " " 
                << col[2] << " ";
        }
    }
    std::vector<Vec3> colors;
    int width, height;
};


int main() {
    Image img(480, 480);
    SobolSampler sampler;
    int cnt = 512;
    for(int i = 0; i < cnt; i++) {
        int y = img.height * sampler.sample(i, 0);
        int x = img.width * sampler.sample(i, 1);
        y = std::min(y, img.height-1);
        x = std::min(x, img.width-1);
        img[img.idx(x,y)]=Vec3(1.0f);
    }
    img.write_ppm("sobol_sample.ppm", false);
}