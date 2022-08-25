/*
c++ cubic_b_spline.cpp -std=c++11 -o cubic_b_spline -lopencv_core -lopencv_highgui -lopencv_imgcodecs -lopencv_imgproc && ./cubic_b_spline
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

float b3(float t) {
    if(t >= 0 && t <= 1)return 1.0/6.0 * std::pow(t, 3);
    if(t >= 1 && t <= 2)return 1.0/6.0 * (-3 * std::pow(t-1, 3) + 3 * std::pow(t-1, 2) + 3 * (t-1) + 1);
    if(t >= 2 && t <= 3)return 1.0/6.0 * ( 3 * std::pow(t-2, 3) - 6 * std::pow(t-2, 2)             + 4);
    if(t >= 3 && t <= 4)return 1.0/6.0 * (-1 * std::pow(t-3, 3) + 3 * std::pow(t-3, 2) - 3 * (t-3) + 1);
    return 0;
}

cv::Point2f gamma(const std::vector<cv::Point2f>& points, float t) {
    cv::Point2f res(0, 0);
    int n = points.size() - 1;
    for(int i = 0; i <= n; i++) {
        res += points[i] * b3(t - i);
    }
    return res;
}

// t >= 3
cv::Point2f foo(const std::vector<cv::Point2f>& points, float t) {
    int j = int(t);
    cv::Mat Gb = (cv::Mat_<float>(2, 4) << 
                  points[j].x, points[j-1].x, points[j-2].x, points[j-3].x,
                  points[j].y, points[j-1].y, points[j-2].y, points[j-3].y);
    cv::Mat MBs = 1.0/6.0 * (cv::Mat_<float>(4, 4) << 
                  0,  0, 0,  1,
                  1,  3, 3, -3,
                  4,  0, -6, 3,
                  1, -3, 3, -1);
    cv::Mat T = (cv::Mat_<float>(4, 1) << 
                  1, (t-j), std::pow(t-j,2), std::pow(t-j,3));
    cv::Mat tmp = Gb * MBs * T;
    cv::Point2f p;
    p.x = tmp.at<float>(0); p.y = tmp.at<float>(1);
    return p;
}

void cubic_b_spline2(const std::vector<cv::Point2f>& points, cv::Mat& window) {
    int n = points.size() - 1;
    float eps = 1e-3;
    for(float t = 3.0; t <= n+1-eps; t += eps) {
        cv::Point2f p = foo(points, t);
        // std::cout << "t, p: " << t << ", " << p << "\n";
        window.at<cv::Vec3b>(p.y, p.x) = cv::Vec3b(255, 0, 0);
    }
}

void cubic_b_spline(const std::vector<cv::Point2f>& points, cv::Mat& window) {
    int n = points.size() - 1;
    float eps = 1e-3;
    for(float t = 3.0; t <= n+1-eps; t += eps) {
        cv::Point2f p = gamma(points, t);
        // std::cout << "t, p: " << t << ", " << p << "\n";
        window.at<cv::Vec3b>(p.y, p.x) = cv::Vec3b(255, 0, 0);
    }
}

int main() 
{
    cv::Mat window = cv::Mat(700, 700, CV_8UC3, cv::Vec3b(255,255,255));
    cv::cvtColor(window, window, cv::COLOR_BGR2RGB);
    cv::namedWindow("Cubic B Spline", cv::WINDOW_AUTOSIZE);

    cv::setMouseCallback("Cubic B Spline", mouse_handler, nullptr);

    int key = -1;
    while (key != 27 && key != 's') 
    {
        for (auto &point : control_points) 
        {
            cv::circle(window, point, 1, {0,0,0}, 3);
        }

        if (key == 'd' && control_points.size() >= 4) 
        {
            cubic_b_spline(control_points, window);
            std::vector<cv::Point2f>().swap(control_points);
        }

        cv::imshow("Cubic B Spline", window);
        key = cv::waitKey(20);
    }
    if(key == 's')cv::imwrite("cubic_b_spline.ppm", window);
    return 0;
}

