/*
clang++ glm_test.cpp -o glm_test && ./glm_test
*/

#include<glm/glm.hpp>
#include<glm/ext.hpp>
#include<iostream>
#include<vector>
using std::vector;
using glm::vec3;
using glm::vec4;
using glm::mat4;

void print(const glm::mat4& mat) {
    for(int i = 0; i < 4; i++) {
        for(int j = 0; j < 4; j++) {
            printf("%7.3f,", mat[i][j]);
        }std::cout << "\n";
    } 
}

void print(const glm::vec4& vec) {
    for(int i = 0; i < 4; i++) {
        printf("%7.3f,", vec[i]);
    }std::cout << "\n";
}

void print(const glm::vec3& vec) {
    for(int i = 0; i < 3; i++) {
        printf("%7.3f,", vec[i]);
    }printf("\n");
}
void print(const glm::i32vec2& vec) {
    for(int i = 0; i < 2; i++) {
        printf("%7d,", vec[i]);
    }printf("\n");
}


vector<vec4> func(vec4 v1, vec4 v2) {
    vector<vec4> res;
    res.push_back(v1);
    res.push_back(v2);
    for(int i = 0; i < 3; i++) {
        float f1, f2;
        for(int j = -1; j <= 1; j+= 2) {
            f1 = res[0].w + j * res[0][i];
            f2 = res[1].w + j * res[1][i];
            if(f1 < 0 && f2 < 0) {
                return vector<vec4>();
            }
            if(f1 >= 0 && f2 >= 0)continue;
            // w+x = 0: a = (w1 + x1)/((w1+x1)-(w2+x2))
            // w-x = 0: a = (-w1+x1)/((-w1+x1)-(-w2+x2))
            float reciprocal = (j * res[0].w + res[0][i] - (j * res[1].w + res[1][i]));
            if(reciprocal == 0)continue;    // 平行
            float alpha = (j * res[0].w + res[0][i]) 
                        / reciprocal;
            // printf("%f\n", alpha);
            vec4 v = (1 - alpha) * res[0] + alpha * res[1];
            if(f1 < 0)res[0] = v;
            else res[1] = v;
            // print(res[0]);
            // print(res[1]);
            // printf("\n");
        }
    }
    return res;
}

vec3 cal_bary3D(const vec3& v, vec3 vecs[]) {
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

int main() {
    float n = 0.3f, f = 100.0f, fov = 70.0f;
    float aspect = 2;
    glm::mat4 projection = glm::perspective(glm::radians(fov), aspect, n, f);
    vec4 v1 = vec4(-2, -2, -3, 1);
    vec4 v2 = vec4( 2, -2,  3, 1);
    // mat4 projection = mat4(n/r,   0,           0,            0,
    //                          0, n/t,           0,            0,
    //                          0,   0,-(f+n)/(f-n), -2*f*n/(f-n),
    //                          0,   0,          -1,            0);
    print(projection);
    printf("\n");
    v1 = projection * v1;
    v2 = projection * v2;
    // vec4 v1 = vec4(0, -1, -2, 1) * projection;
    // vec4 v2 = vec4(0, -1,  2, 1) * projection;
    print(v1);
    print(v2);
    printf("\n");

    vector<vec4> res = func(v1, v2);
    printf("size: %d\n", (int)res.size());
    if(res.size() == 2) {
        print(res[0]);
        print(res[1]);
        // print(res[0]/res[0].w);
        // print(res[1]/res[1].w);
    }

    vec3 va(0, 1, 0);
    vec3 vb(-1, 0, 0);
    vec3 vc(1, 0, 0);
    vec3 vp(0, 2, 0);
    vec3 vecs[3] = {va, vb, vc};
    vec3 weights = cal_bary3D(vp, vecs);
    print(weights);
    return 0;
}
