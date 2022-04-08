/*
clang++ gamma.cpp -o gamma && ./gamma
*/
#include<iostream>
#include<math.h>
#include<stdlib.h>
#define W 1024
#define H 512

int gamma_encode(float radiance, bool encode) {
    if(encode)
        return (int)(pow(std::min(1.0f, std::max(0.0f, radiance)), 1.0f/2.2f) * 255.0f);
    else
        return (int)(radiance * 255.0f);
}

void write(const char* filename) {
    FILE* fp = fopen(filename, "w");
    fprintf(fp, "P3 %d %d 255\n", W, H);
    for(int i = 0; i < H; i++) {
        float r = (float)(i/2) / 255.0f;
        for(int j = 0; j < W; j++) {
            int c1 = gamma_encode(r, true);
            int c2 = gamma_encode(r, false);
            if(j < W/2) {
                fprintf(fp, "%d %d %d ", c1, c1, c1);
            }
            else
                fprintf(fp, "%d %d %d ", c2, c2, c2);
        }
    }
    fclose(fp);
    return;
}

int main() {
    const char* filename = "gamma_vs_nogamma.ppm";
    write(filename);
    return 0;
}