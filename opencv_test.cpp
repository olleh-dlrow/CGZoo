/*
clang++ opencv_test.cpp -std=c++11 -o opencv_test -L/usr/local/lib/opencv4.5.5 -lopencv_core -lopencv_highgui -lopencv_imgcodecs && ./opencv_test
*/
#include<stdio.h>
#include<opencv2/opencv.hpp>

int main() {
    const char* path = "image.ppm";
    cv::Mat img = cv::imread(path);
    int key;
    while(key != 'q') {
        cv::imshow("origin", img);
        key = cv::waitKey(10);
    }
    return 0;
}