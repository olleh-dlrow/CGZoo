/*
c++ path_tracing.cpp -O3 -std=c++20 -o path_tracing && ./path_tracing
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

/////////////////////////////
// DECLARATION
/////////////////////////////
struct Mesh;
struct Material;
struct Light;
using Float = float;

/////////////////////////////
// CONST VARIABLE
/////////////////////////////
const Float PI = std::acos(-1);
const Float EPSILON = 1e-5;
const Float FLOAT_INF = std::numeric_limits<Float>::infinity();
const int INT32_INF = std::numeric_limits<int>::infinity();
/////////////////////////////
// UTILS
/////////////////////////////

Float get_uniform01()  {
    static std::uniform_real_distribution<Float> dist;
    static std::random_device dev;
    static std::mt19937 rng(dev());
    return dist(rng);
}

Float deg2rad(const Float& deg) { return deg * PI / 180.0; }

int gamma_encode(Float radiance, bool gamma) {
    if(gamma)
        return (int)(pow(std::min(1.0f, std::max(0.0f, radiance)), 1.0f/2.2f) * 255.0f);
    else 
        return (int)(radiance * 255.0f);
}

Float inv_sqrt (Float x)
{
    Float xhalf = 0.5f*x;
    int i = *(int*)&x;
    i = 0x5f3759df - (i >> 1); 
    x = *(Float*)&i;
    x = x*(1.5f - xhalf*x*x); 
    return x;
}

template<typename T>
T clamp(T v, T minv, T maxv) {
    return std::max(std::min(v, maxv), minv);
}

/////////////////////////////
// BASE
/////////////////////////////

struct Color24 {
    unsigned char r = 0;
    unsigned char g = 0;
    unsigned char b = 0;

    Color24 operator+(const Color24& rhs) const {
        return Color24{
            (unsigned char)clamp((int)r + rhs.r, 0, 255),
            (unsigned char)clamp((int)g + rhs.g, 0, 255),
            (unsigned char)clamp((int)b + rhs.b, 0, 255)
        };
    };
};

template<>
Color24 clamp(Color24 v, Color24 min, Color24 max) {
    return Color24{
        clamp(v.r, min.r, max.r),
        clamp(v.g, min.g, max.g),
        clamp(v.b, min.b, max.b),
    };
}

struct Vec2 {
    Float x = 0;
    Float y = 0;

    Float operator[](int index) {
        switch(index) {
            case 0:
                return x;
            case 1:
                return y;
            default:
                assert(0);
        }
        return FLOAT_INF;
    }

    Vec2 operator+(const Vec2& rhs) const {
        return {x + rhs.x, y + rhs.y};
    }
};

Vec2 operator*(Float lhs, const Vec2& rhs) {
    return {lhs * rhs.x, lhs * rhs.y};
}

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

    Color24 to_color() const {
        Vec3 v = (*this + Vec3(1.0f)) * 0.5f * 255;
        return Color24{
            (unsigned char)clamp(v.x, 0.0f, 255.0f),
            (unsigned char)clamp(v.y, 0.0f, 255.0f),
            (unsigned char)clamp(v.z, 0.0f, 255.0f),
        };
    }

    Vec3 normalize() const {
        return *this * inv_magnitude();
    };

    Float inv_magnitude() const {
        return inv_sqrt(x*x + y*y + z*z);
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
        return FLOAT_INF;
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

struct Ray {
    Ray(const Vec3& o, const Vec3& d):
        org(o), dir(d) {
        dir_inv = 1. / dir;
    }

    Vec3 value(Float t) const {
        return org + t * dir;
    }

    Vec3 org;
    Vec3 dir;
    Vec3 dir_inv;
};

struct Transfrom {
    Transfrom() {

    }
    Transfrom(const Vec3& p, const Vec3& s):
        position(p), scale(s){

    }
    Vec3 position = Vec3(0.0f);
    Vec3 scale = Vec3(1.0f);
};

struct Bound {
    Bound():min(0, 0, 0), max(1, 1, 1) {

    }

    Bound(Vec3 v1, Vec3 v2) {
        min = Vec3::min(v1, v2);
        max = Vec3::max(v1, v2);
    }

    bool is_intersect(Ray& ray) const {
        Float t_enter = -FLOAT_INF;
        Float t_exit = FLOAT_INF;
        for(int i = 0; i < 3; i++) {
            if(std::fabs(ray.dir[i]) == FLOAT_INF) {
                continue;
            }
            Float t_min = (min[i] - ray.org[i]) * ray.dir_inv[i];
            Float t_max = (max[i] - ray.org[i]) * ray.dir_inv[i];

            if(ray.dir[i] < 0) {
                std::swap(t_min, t_max);
            }
            
            t_enter = std::max(t_enter, t_min);
            t_exit = std::min(t_exit, t_max);
        }
        return t_enter <= t_exit && t_exit >= 0;
    }

    std::string to_string() const {
        return "[" + min.to_string() + ", " + max.to_string() + "]";
    }

    Vec3 centroid() const {
        return (min + max) / 2.0f;
    }

    Bound union_bound(const Bound& rhs) const {
        return Bound(Vec3::min(min, rhs.min) - Vec3(EPSILON), 
                     Vec3::max(max, rhs.max) + Vec3(EPSILON)
                    );
    }

    Bound union_vec3(const Vec3& rhs) const {
        return Bound(
            Vec3::min(min, rhs) - Vec3(EPSILON),
            Vec3::max(max, rhs) + Vec3(EPSILON)
        );
    }

    int max_axis() const {
        int indices[3] = {0,1,2};
        Float axises[3] = {max.x - min.x, max.y - min.y, max.z - min.z};
        std::sort(indices, indices + 3, [&axises](int i, int j){
            return axises[i] > axises[j];
        });
        return indices[0];
    }


    Vec3 min;
    Vec3 max;
};

struct Intersection {
    bool happened = false;
    Vec3 p = Vec3();
    Vec3 normal = Vec3();
    Float t = -FLOAT_INF;
    Material* mat = nullptr;
    Light* light = nullptr;
};

struct ShapeSample {
    Intersection inter;
    Float pdf;
};

/////////////////////////////
// SHAPE
/////////////////////////////

struct Triangle {
    Triangle(Mesh* mesh, int _face_idx);

    Float get_area() const {
        return area;
    }

    Bound get_bound() const {
        return Bound(vertices[0], vertices[1]).union_vec3(vertices[2]);
    }

    Vec3 get_normal() const {
        return normal;
    }

    const Vec3* get_vertices() const {
        return vertices;
    }

    ShapeSample sample(Vec2 u) const {
        if(u.x + u.y > 1) {
            u.x = 1 - u.x;
            u.y = 1 - u.y;
        }
        Float b[3] = {u.x, u.y, 1 - u.x - u.y};
        Intersection inter;
        inter.happened = true;
        inter.p = vertices[0] * b[0] + vertices[1] * b[1] + vertices[2] * b[2];
        inter.normal = normal;
        ShapeSample s {
            inter, 1.0f / area
        };
        return s;
    }

    Intersection intersect(Ray& ray) const {
        Intersection inter;
        if(ray.dir.dot(normal) > 0) {
            return inter;
        }
        Float u, v, t = 0;
        Vec3 s1 = ray.dir.cross(e2);//pvec
        Float det = s1.dot(e1);
        if(fabs(det) < EPSILON) {
            return inter;
        }
        Float det_inv = 1. / det;
        Vec3 s = ray.org - vertices[0];//tvec
        Vec3 s2 = s.cross(e1);//qvec
        u = s1.dot(s) * det_inv;
        if(u < 0 || u > 1) {
            return inter;
        }
        v = s2.dot(ray.dir) * det_inv;
        if(v < 0 || u + v > 1) {
            return inter;
        }
        t = s2.dot(e2) * det_inv;
        if(t <= 0.) {
            return inter;
        }
        inter.happened = true;
        inter.normal = normal;
        inter.t = t;
        inter.p = ray.value(t);

	    return inter;
    }

    // truly intersect with this triangle
    bool intersect(Ray& ray, Intersection& inter) const {
        if(ray.dir.dot(normal) > 0) {
            return false;
        }
        Float u, v, t = 0;
        Vec3 s1 = ray.dir.cross(e2);//pvec
        Float det = s1.dot(e1);
        if(fabs(det) < EPSILON) {
            return false;
        }
        Float det_inv = 1. / det;
        Vec3 s = ray.org - vertices[0];//tvec
        Vec3 s2 = s.cross(e1);//qvec
        u = s1.dot(s) * det_inv;
        if(u < 0 || u > 1) {
            return false;
        }
        v = s2.dot(ray.dir) * det_inv;
        if(v < 0 || u + v > 1) {
            return false;
        }
        t = s2.dot(e2) * det_inv;
        if(t <= 0.) {
            return false;
        }

        if(!inter.happened || t < inter.t) {
            inter.happened = true;
            inter.normal = normal;
            inter.t = t;
            inter.p = ray.value(t);
            return true;
        }
        return false;
    }


    Float area;
    Vec3 normal;
    Vec3 vertices[3];
    Vec3 e1, e2;
    int face_index;
};

/////////////////////////////
// MATERIAL
/////////////////////////////

enum MatType {
    DIFFUSE, 
    SPECULAR
};

struct Material {
    Material(
        MatType t = DIFFUSE,
        const Vec3& e = Vec3(),
        const Vec3& kd = Vec3(),
        const Vec3& ks = Vec3()
    ):
        type(t), emission(e), Kd(kd), Ks(ks) {

    }

    Vec3 get_bsdf(const Vec3& wi, const Vec3& wo, const Vec3& N) const {
        switch(type) {
            case DIFFUSE:
                return Kd / PI;
            case SPECULAR:
                return Vec3();
        }
        return Vec3();
    }

    MatType type;
    Vec3 emission;
    Vec3 Kd, Ks;
};

/////////////////////////////
// MESH
/////////////////////////////

struct Mesh {
    Mesh(const std::string& filename, const Transfrom& trans=Transfrom()) {
        load(filename, trans);
    }

    Bound get_bound() const {
        return bound;
    }

    void load(const std::string& filename, const Transfrom& trans) {
        std::ifstream file(filename);
        if(!file) {
            std::cout << "can't load file: " << filename << std::endl;
            return;
        }
        
        Vec3 min_vert = Vec3(FLOAT_INF);
        Vec3 max_vert = Vec3(-FLOAT_INF);
        int cur_face_index = 0;
        std::string line;
        std::string delim(" \r\n");
        while(std::getline(file, line)) {
            if(line[0] == 'v') {
                int s = 1, e = 0;
                s = line.find_first_not_of(delim, s);
                e = line.find_first_of(delim, s);
                std::array<Float, 3> vertex;
                int i = 0;
                while(s != std::string::npos) {
                    vertex[i++] = std::stof(line.substr(s, e - s));
                    s = line.find_first_not_of(delim, e);
                    e = line.find_first_of(delim, s);
                }
                Vec3 vert = Vec3{vertex[0], vertex[1], vertex[2]} * trans.scale
                            + trans.position;
                min_vert = Vec3::min(min_vert, vert); 
                max_vert = Vec3::max(max_vert, vert);              
                vertices.emplace_back(vert);
            } else if(line[0] == 'f') {
                int s = 1, e = 0;
                s = line.find_first_not_of(delim, s);
                e = line.find_first_of(delim, s);
                while(s != std::string::npos) {
                    indices.emplace_back(std::stoi(line.substr(s, e - s)) - 1);
                    s = line.find_first_not_of(delim, e);
                    e = line.find_first_of(delim, s);                    
                }
                tris.emplace_back(Triangle(this, cur_face_index++));
            }
        }
        bound = Bound(min_vert, max_vert);
        return;
    }    

    std::string to_string() const {
        std::stringstream ss;
        ss << "vertices count: " << vertices.size() << std::endl;
        ss << "indices count: " << indices.size() << std::endl;
        ss << "triangles count: " << tris.size() << std::endl;
        return ss.str();
    }

    Mesh(Mesh&& val):
        vertices(std::move(val.vertices)),
        normals(std::move(val.normals)),
        uvs(std::move(val.uvs)),
        indices(std::move(val.indices)),
        tris(std::move(val.tris)),
        mat(val.mat) {
            val.mat = nullptr;
    }

    Mesh& operator=(Mesh&& rhs) {
        Mesh(std::move(rhs)).swap(*this);
        return *this;
    }

    void swap(Mesh& rhs) {
        vertices.swap(rhs.vertices);
        normals.swap(rhs.normals);
        uvs.swap(rhs.uvs);
        indices.swap(rhs.indices);
        tris.swap(rhs.tris);
        std::swap(mat, rhs.mat);
    }

    Bound bound;
    std::vector<Vec3> vertices;
    std::vector<Vec3> normals;
    std::vector<Vec2> uvs;
    std::vector<int> indices;
    std::vector<Triangle> tris;

    Material* mat = nullptr;
};

/////////////////////////////
// LIGHT
/////////////////////////////

struct LightSample {
    Intersection inter;
    Float pdf;
    Vec3 L;
};

struct Light {
    Light(Triangle* _tri, Material* _mat):
        tri(_tri), mat(_mat) {
        
    }

    Vec3 get_L() const {
        return mat->emission;
    }

    LightSample sample(Vec2 u) const {
        ShapeSample s = tri->sample(u);
        LightSample ls;
        s.inter.mat = mat;

        ls.inter = s.inter;
        ls.pdf = s.pdf;
        ls.L = get_L();
        return ls;
    }

    Triangle* tri;
    Material* mat;
};

/////////////////////////////
// PRIMITIVE
/////////////////////////////

struct Primitive {
    Primitive(Triangle* _tri, Material* _mat, Light* li=nullptr):
        tri(_tri), mat(_mat), light(li) {

    }

    Bound get_bound() const {
        return tri->get_bound();
    }

    Intersection intersect(Ray& ray) {
        Intersection inter = tri->intersect(ray);
        if(inter.happened) {
            inter.mat = mat;
            inter.light = light;
        }
        return inter;
    }

    void intersect(Ray& ray, Intersection& inter) {
        bool is_inter = tri->intersect(ray, inter);
        if(is_inter) {
            inter.mat = mat;
            inter.light = light;
        }
    }

    ShapeSample sample(Vec2 u) const {
        ShapeSample s = tri->sample(u);
        s.inter.mat = mat;
        return s;
    }

    Triangle* tri;
    Material* mat;
    Light* light = nullptr;
};

/////////////////////////////
// BVH
/////////////////////////////

struct BVHNode {
    BVHNode(const Bound& b):bound(b) {

    }

    bool is_bound_intersect(Ray& ray) {
        return bound.is_intersect(ray);
    }

    void inner_intersect(Ray& ray, Intersection& inter) {
        for(auto&& prim : prims) {
            prim->intersect(ray, inter);
        }
    }

    bool is_leaf() const {
        return !left && !right;
    }

    std::string to_string() const {
        return " bound: " + bound.to_string() 
                + " prim size: " + std::to_string(prims.size());
    }

    Bound bound;
    BVHNode* left = nullptr, *right = nullptr;
    std::vector<Primitive*> prims;
};

struct BVH {
    BVH(int _max_prims_in_node=1, std::vector<Primitive>&& _prims=std::vector<Primitive>())
        :max_prims_in_node(_max_prims_in_node), prims(std::move(_prims)) {
        
    }

    void intersect(Ray& ray, Intersection& inter) {
        recursive_intersect(root, ray, inter);
    }

    void recursive_intersect(BVHNode* node, Ray& ray, Intersection& inter) {
        if(!node->is_bound_intersect(ray)) {
            return;
        } 
        if(node->is_leaf()) {
            node->inner_intersect(ray, inter);
        } else {
            recursive_intersect(node->left, ray, inter);
            recursive_intersect(node->right, ray, inter);
        }
    }
    
    void build() {
        std::vector<Primitive*> ptrs(prims.size());
        for(int i = 0; i < prims.size(); i++) {
            ptrs[i] = &prims[i];
        }
        root = recursive_build(std::move(ptrs));
    }

    BVHNode* recursive_build(std::vector<Primitive*> ptrs) {
        // calculate bound
        Bound b = ptrs[0]->get_bound();
        for(int i = 1; i < ptrs.size(); i++) {
            b = b.union_bound(ptrs[i]->get_bound());
        }
        BVHNode* node = new BVHNode(b);
        if(ptrs.size() > max_prims_in_node) {
            // split
            int dim = b.max_axis();
            std::sort(ptrs.begin(), ptrs.end(), [&dim](Primitive* p1, Primitive* p2){
                return p1->get_bound().centroid()[dim] < p2->get_bound().centroid()[dim];
            });            
            auto mid = ptrs.begin() + ptrs.size() / 2;
            std::vector<Primitive*> left_split(ptrs.begin(), mid);
            std::vector<Primitive*> right_split(mid, ptrs.end());
            node->left = recursive_build(std::move(left_split));
            node->right = recursive_build(std::move(right_split));
        } else {
            node->prims = std::move(ptrs);
        }
        return node;
    }

    std::string to_string() const {
        std::stringstream ss;
        std::queue<BVHNode*> qu;
        qu.push(root);
        int idx = 0;
        while(!qu.empty()) {
            auto cur = qu.front();
            qu.pop();
            ss << idx++ << cur->to_string() << std::endl;
            if(!cur->is_leaf()) {
                qu.push(cur->left);
                qu.push(cur->right);
            }
        }
        return ss.str();
    }

    ~BVH() {
        // implementation
    }

    int max_prims_in_node;
    BVHNode* root;
    std::vector<Primitive> prims;
};

/////////////////////////////
// SCENE
/////////////////////////////

struct Scene {
    void intersect(Ray& ray, Intersection& inter) {
        bvh.intersect(ray, inter);
    }

    void add_mesh(Mesh* mesh) {
        meshes.emplace_back(mesh);
    }

    void add_light(Light* light) {
        lights.emplace_back(light);
    }

    void build_bvh() {
        std::vector<Primitive> prims;
        for(int i = 0; i < meshes.size(); i++) {
            auto& mesh = meshes[i];
            bool has_light = mesh2light.count(i) > 0;
            for(auto&& tri : mesh->tris) {
                prims.emplace_back(
                    &tri,
                    mesh->mat
                );
                if(has_light) {
                    prims.back().light = mesh2light[i];
                }
            }
        }
        bvh = BVH(1, std::move(prims));
        bvh.build();
        return;
    }

    std::vector<Mesh*> meshes;
    std::vector<Light*> lights;
    std::unordered_map<int, Light*> mesh2light;

    BVH bvh;
};

/////////////////////////////
// SAMPLER
/////////////////////////////

struct Sampler {
    Sampler() {

    }

    void update(unsigned spp_index) {
        this->spp_index = spp_index;
    }

    virtual Float sample(unsigned dim) const = 0;

    unsigned spp_index;
};

struct SobolSampler: public Sampler {
    SobolSampler() {

    }

    Float sample(unsigned dim) const override {
        return P[spp_index][dim];
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
    const static unsigned MAXN = 2048;
    const static unsigned MAXD = 7;
};

double** SobolSampler::P = SobolSampler::sobol_points(
    SobolSampler::MAXN, SobolSampler::MAXD, "new-joe-kuo-6.21201");

/////////////////////////////
// CAMERA
/////////////////////////////

struct Camera {
    Camera(Vec3 e, int w, int h, Float _fov, Float _near=1.0f):
        eye(e), width(w), height(h), aspect(w/(Float)h), fov(_fov), near(_near){
       scale = std::tan(deg2rad(fov) * 0.5f);
    }

    Ray cast_ray(const Vec2& pixel, Sampler& sampler) {
        Float x = 
        (2 * (pixel.x + sampler.sample(0)) / (Float)width - 1.0f) 
        * aspect * scale;
        Float y = 
        (2 * (pixel.y + sampler.sample(1)) / (Float)height - 1.0f)
        * scale;
        Vec3 dir(x, y, near);
        return Ray(eye, dir.normalize());
    }

    Float scale;
    Vec3 eye;
    Float aspect;
    Float fov;
    Float near = 1.0f;
    int width, height;
};

/////////////////////////////
// INTEGRATOR
/////////////////////////////

struct Integrator {
    Integrator(Scene* s, Float rr=0.8):
        scene(s), russian_roulette(rr) {

    }

    Vec3 Li(Ray& ray, Sampler& sampler) {
        bool sampled_light = false;
        return recursive_Li(ray, sampler, sampled_light);
    }

    Vec3 recursive_Li(Ray& ray, Sampler& sampler, bool& sampled_light) {
        Intersection inter;
        scene->intersect(ray, inter);
        if(!inter.happened) {
            sampled_light = false;
            return Vec3(0);
        }

        // intersect with light
        if(inter.light) {
            sampled_light = true;
            return inter.light->get_L();
        }

        Vec3 wo = -1 * ray.dir.normalize();
        // give some offset to avoid intersection itself
        Float eps = 0.001f;
        Vec3 offset_inter_p = inter.p + eps * inter.normal;

        // contribution from light source
        Vec3 L_dir = Vec3(0);
        for(auto&& li : scene->lights) {
            // check block
            LightSample s = li->sample(
                {get_uniform01(), 
                 get_uniform01()}
            );

            Vec3 vp2li = s.inter.p - offset_inter_p;
            Ray p2li(offset_inter_p, vp2li.normalize());
            Intersection block;
            scene->intersect(p2li, block);
            // if blocked, continue
            if((block.p - s.inter.p).sqr_magnitude() > eps) {
                continue;
            }
            Vec3 L_i = s.L;
            Vec3 wi = p2li.dir.normalize();
            Vec3 f_r = inter.mat->get_bsdf(wi, wo, inter.normal);
            Float cos_theta = wi.dot(inter.normal);
            Float cos_alpha = s.inter.normal.dot(-1 * wi);
            Float pdf_light = s.pdf;
            Float sqr_dist = vp2li.sqr_magnitude();
            L_dir += L_i * f_r * cos_alpha * cos_theta / sqr_dist / pdf_light;
        }

        // contribution from other reflections
        Vec3 L_indir = Vec3(0);
        Float r = get_uniform01();
        if(r < russian_roulette) {
            Vec3 wi = sample_uniform_hemisphere({
                get_uniform01(),
                get_uniform01()
            });
            if(wi.dot(inter.normal) * wo.dot(inter.normal) < 0) {
                wi = -1 * wi;
            }
            Ray next(offset_inter_p, wi);
            L_indir = recursive_Li(next, sampler, sampled_light);
            if(sampled_light) {
                L_indir = Vec3(0.0f);
            } else {
                Float pdf_hemi = uniform_hemisphere_pdf();
                Float cos_theta = wi.dot(inter.normal);
                Vec3 f_r = inter.mat->get_bsdf(wi, wo, inter.normal);
                L_indir = L_indir * f_r * cos_theta / pdf_hemi / russian_roulette;
            }
        }
        sampled_light = false;
        return L_dir + L_indir;
    }

    static Float uniform_hemisphere_pdf() {
        // area: 1 / (2pi)
        return 0.5f / PI;
    }

    static Vec3 sample_uniform_hemisphere(Vec2 u) {
        Float z = u[0];
        Float r = std::sqrt(1 - z*z);
        Float phi = 2 * PI * u[1];
        return {r * std::cos(phi), r * std::sin(phi), z};        
    }

    Scene* scene;
    Float russian_roulette;
};

/////////////////////////////
// IMAGE
/////////////////////////////

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

/////////////////////////////
// RENDERER
/////////////////////////////

struct Renderer {
    Renderer(Integrator* inte, Camera* cam, Image* img):
        integrator(inte), camera(cam), image(img) {
        
    }

    void render(int spp) {

        std::cout << "start render...\n";
        std::cout << "spp: " << spp << "\n";
        std::cout << "rr: " << integrator->russian_roulette << "\n";
        std::cout << "threads: \n";
        // work threads
        int nums_threads = std::thread::hardware_concurrency();
        int num_lines = image->height / nums_threads + 1;
        std::vector<std::thread> work_threads;
        std::vector<long> finish_cnts(nums_threads, 0);

        auto start = std::chrono::system_clock::now();
        for(int i = 0; i < nums_threads; i++) {
            int y0 = i * num_lines;
            int y1 = std::min(y0 + num_lines, image->height);
            std::cout << "thread " << i << ": " << y0 << " => " 
                << y1 << std::endl;
            work_threads.emplace_back(std::thread(
                [this, &spp, &finish_cnts]
                (int sy, int ey, int thread_idx){
                SobolSampler sampler;
                for(int y = sy; y < ey; y++) {
                    for(int x = 0; x < image->width; x++) {
                        Vec2 pixel{(float)x, (float)y};
                        for(int i = 0; i < spp; i++) {
                            sampler.update(i);
                            Ray ray = camera->cast_ray(pixel, sampler);
                            int index = pixel.y * image->width + pixel.x;
                            Vec3 L = integrator->Li(ray, sampler) / spp;
                            (*image)[image->idx(-x, y)] += L;
                            finish_cnts[thread_idx]++;
                        }
                    }
                }
            }, y0, y1, i));
        }

   
            char roll_sign[3] = {'-','\\','/'};
            int cur_roll_idx = 0;
            long sum = image->width * image->height * spp;
            
        for(;;) {
            long cur_sum = 0;
            for(auto&& cnt : finish_cnts) {
                cur_sum += cnt;
            }
            auto cur_time = std::chrono::system_clock::now();
            long seconds = std::chrono::duration_cast<std::chrono::seconds>(cur_time - start).count();
            
            std::printf("\r%.2f%% time: %ld:%02ld:%02ld %c",
                    (Float)cur_sum / sum * 100, seconds / 3600, (seconds / 60) % 60, seconds % 60, roll_sign[cur_roll_idx]);
            std::fflush(stdout);
            cur_roll_idx = (cur_roll_idx + 1) % 3;                
            if(cur_sum >= sum) {
                break;
            } 
            using namespace std::chrono_literals;
            std::this_thread::sleep_for(1s);
        }

        for (auto &&t : work_threads)
        {
            t.join();
        }
        std::cout << "\nrender finish\n";
    }

    Integrator* integrator;
    Camera* camera;
    Image* image;
};

/////////////////////////////
// IMPLEMENTATION
/////////////////////////////

Triangle::Triangle(Mesh* mesh, int _face_idx): 
    face_index(_face_idx) {
    int *idx = &mesh->indices[3 * face_index];
    vertices[0] = mesh->vertices[idx[0]]; 
    vertices[1] = mesh->vertices[idx[1]];
    vertices[2] = mesh->vertices[idx[2]];
    e1 = vertices[1] - vertices[0];
    e2 = vertices[2] - vertices[0];
    normal = e1.cross(e2).normalize();
    // area
    area = e1.cross(e2).magnitude() * 0.5f;
}

/////////////////////////////
// MAIN
/////////////////////////////


int main() {
    Scene scene;
    // load materials
    std::unordered_map<std::string, std::unique_ptr<Material>> materials;
    materials.insert({"Red",   std::move(std::make_unique<Material>(DIFFUSE, Vec3(), Vec3(0.63f, 0.065f, 0.05f)))});
    materials.insert({"Green", std::move(std::make_unique<Material>(DIFFUSE, Vec3(), Vec3(0.14f, 0.45f, 0.091f)))});
    materials.insert({"White", std::move(std::make_unique<Material>(DIFFUSE, Vec3(), Vec3(0.725f, 0.71f, 0.68f)))});
    materials.insert({"Light", std::move(std::make_unique<Material>(
            DIFFUSE, 
            ( 8.0f * Vec3(0.747f+0.058f, 0.747f+0.258f, 0.747f) 
            + 15.6f * Vec3(0.740f+0.287f,0.740f+0.160f,0.740f) 
            + 18.4f * Vec3(0.737f+0.642f,0.737f+0.159f,0.737f)),
             Vec3(0.65f, 0.65f, 0.65f)
    ))});
    
    // load meshes
    std::vector<std::unique_ptr<Mesh>> meshes;
    std::string bunny = "./models/bunny/bunny.obj";
    std::string prefix = "./models/cornellbox/";
    std::vector<std::string> mesh_files = {
        "floor.obj",
        "left.obj",
        "right.obj",
        "shortbox.obj",
        "tallbox.obj",
        // "../bunny/bunny.obj",
    };
    std::vector<Transfrom> mesh_trans = {
        Transfrom(),
        Transfrom(),
        Transfrom(),
        Transfrom(),
        Transfrom(),
        Transfrom(Vec3(300, 300, 300), Vec3(3000.0f))
    };
    std::vector<std::string> mesh_mats = {
        "White",
        "Red",
        "Green",
        "White",
        "White",
        "White"
    };
    for(int i = 0; i < mesh_files.size(); i++) {
        auto& filename = mesh_files[i];
        meshes.emplace_back(std::make_unique<Mesh>(prefix + filename, mesh_trans[i]));
        meshes[i]->mat = materials[mesh_mats[i]].get();
        scene.add_mesh(meshes[i].get());
    }
    // load lights
    std::vector<std::unique_ptr<Light>> lights;
    std::vector<std::string> light_files = {
        "light.obj",
    };
    std::vector<std::string> light_mats = {
        "Light",
    };
    for(int i = 0, base = meshes.size(); i < light_files.size(); i++) {
        auto& filename = light_files[i];
        meshes.emplace_back(std::make_unique<Mesh>(prefix + filename));
        auto mesh = meshes.back().get();
        auto mat = materials[light_mats[i]].get();
        meshes[base + i]->mat = mat;

        scene.add_mesh(mesh);

        for(auto&& tri : mesh->tris) {
            lights.emplace_back(std::make_unique<Light>(&tri, mat));
            auto light = lights.back().get();
            scene.mesh2light[base + i] = light;

            scene.add_light(light);
        }
    }

    // build bvh
    scene.build_bvh();
    // image
    Image img(784, 784);
    // camera
    Vec3 eye = Vec3(278, 273, -800);
    Camera camera(eye, img.width, img.height, 40);
    // integrator
    Integrator integrator(&scene, 0.8f);
    Renderer renderer(&integrator, &camera, &img);
    int spp = 1;
    renderer.render(spp);
    img.write_ppm("path_tracing" + std::to_string(spp) + "spp.ppm", true);
    return 0;
}
