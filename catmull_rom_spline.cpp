/*
c++ catmull_rom_spline.cpp -std=c++11 -o catmull_rom_spline -lopencv_core -lopencv_highgui -lopencv_imgcodecs -lopencv_imgproc && ./catmull_rom_spline
*/
#include <iostream>
#include <opencv2/opencv.hpp>

std::vector<cv::Point2f> control_points;

void mouse_handler(int event, int x, int y, int flags, void *userdata) 
{
    if (event == cv::EVENT_LBUTTONDOWN) 
    {
        std::cout << "Left button of the mouse is clicked - position (" << x << ", "
        << y << ")" << '\n';
        control_points.emplace_back(x, y);
    }     
}

float p1(float t) { return 0.5 * (-1 * t*t*t + 2 * t*t - 1 * t); }
float p2(float t) { return 0.5 * ( 3 * t*t*t - 5 * t*t + 2 * 1); }
float p3(float t) { return 0.5 * (-3 * t*t*t + 4 * t*t + 1 * t); }
float p4(float t) { return 0.5 * ( 1 * t*t*t - 1 * t*t); }

float bCR(float t) {
    if(t >= -2 && t <= -1)return p4(t + 2);
    if(t >= -1 && t <=  0)return p3(t + 1);
    if(t >=  0 && t <=  1)return p2(t);
    if(t >=  1 && t <=  2)return p1(t - 1);
    return 0;
}

cv::Point2f gamma(const std::vector<cv::Point2f>& points, float t) {
    cv::Point2f res;
    int n = points.size() - 1;
    for(int i = 0; i <= n; i++) {
        res += points[i] * bCR(t - i);
    }
    return res;
}

void catmull_rom_spline(const std::vector<cv::Point2f>& points, cv::Mat& window) {
    int n = points.size() - 1;
    float eps = 1e-3;
    for(float t = 0.0; t <= n; t += eps) {
        cv::Point2f p = gamma(points, t);
        window.at<cv::Vec3b>(p.y, p.x) = cv::Vec3b(255, 0, 0);
    }
}

int main() 
{
    cv::Mat window = cv::Mat(700, 700, CV_8UC3, cv::Vec3b(255,255,255));
    cv::cvtColor(window, window, cv::COLOR_BGR2RGB);
    cv::namedWindow("Catmull-Rom Spline", cv::WINDOW_AUTOSIZE);

    cv::setMouseCallback("Catmull-Rom Spline", mouse_handler, nullptr);

    int key = -1;
    while (key != 27 && key != 's') 
    {
        for (auto &point : control_points) 
        {
            cv::circle(window, point, 1, {0,0,0}, 3);
        }

        if (key == 'd') 
        {
            catmull_rom_spline(control_points, window);
            std::vector<cv::Point2f>().swap(control_points);
        }

        cv::imshow("Catmull-Rom Spline", window);
        key = cv::waitKey(20);
    }
    if(key == 's')cv::imwrite("catmull_rom_spline.ppm", window);
    return 0;
}
