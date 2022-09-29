#include<iostream>
#include<vector>
#include<random>
using namespace std;

float get_uniform01()  {
    static std::uniform_real_distribution<float> dist;
    static std::random_device dev;
    static std::mt19937 rng(dev());
    return dist(rng);
}

int main() {
    int sum = 500000;
    int cnt = sum;
    float aver = 0.0f;
    float D = 0.0f;
    while(cnt--) {
        float v = get_uniform01();
        aver += v;
        D += (v - 0.5f)*(v-0.5f);
    }
    cout << endl << aver / sum << endl << D / sum << endl;
    return 0;
}