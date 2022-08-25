/*
c++ bezier_curve.cpp -std=c++11 -o bezier_curve -lopencv_core -lopencv_highgui -lopencv_imgcodecs -lopencv_imgproc && ./bezier_curve
*/
#include <iostream>
#include <opencv2/opencv.hpp>

std::vector<cv::Point2f> control_points;

void mouse_handler(int event, int x, int y, int flags, void *userdata) 
{
    if (event == cv::EVENT_LBUTTONDOWN && control_points.size() < 4) 
    {
        std::cout << "Left button of the mouse is clicked - position (" << x << ", "
        << y << ")" << '\n';
        control_points.emplace_back(x, y);
    }     
}

void naive_bezier(const std::vector<cv::Point2f> &points, cv::Mat &window) 
{
    auto &p_0 = points[0];
    auto &p_1 = points[1];
    auto &p_2 = points[2];
    auto &p_3 = points[3];

    for (double t = 0.0; t <= 1.0; t += 0.001) 
    {
        auto point = std::pow(1 - t, 3) * p_0 + 3 * t * std::pow(1 - t, 2) * p_1 +
                 3 * std::pow(t, 2) * (1 - t) * p_2 + std::pow(t, 3) * p_3;

        window.at<cv::Vec3b>(point.y, point.x)[2] = 0;
    }
}

cv::Point2f lerp_v2f(const cv::Point2f& a, const cv::Point2f& b, float t)
{
    return a + (b - a) * t;
}


cv::Point2f recursive_bezier(const std::vector<cv::Point2f> &control_points, float t) 
{
    if (control_points.size() == 1)
    {
        return control_points[0];
    }


    std::vector<cv::Point2f> lerp_points;
    for (size_t i = 1;i<control_points.size();i++)
    {
        lerp_points.push_back(lerp_v2f(control_points[i - 1], control_points[i], t));
    }
    return recursive_bezier(lerp_points,t);
}

void bezier(const std::vector<cv::Point2f> &control_points, cv::Mat &window) 
{
    // TODO: Iterate through all t = 0 to t = 1 with small steps, and call de Casteljau's 
    // recursive Bezier algorithm.
    for (double t = 0.0;t<=1.0;t+=0.001)
    {
        auto point = recursive_bezier(control_points, t);
        window.at<cv::Vec3b>(point.y, point.x)[1] = 0;
    }

}

int main() 
{
    cv::Mat window = cv::Mat(700, 700, CV_8UC3, cv::Vec3b(255,255,255));
    cv::cvtColor(window, window, cv::COLOR_BGR2RGB);
    cv::namedWindow("Bezier Curve", cv::WINDOW_AUTOSIZE);

    cv::setMouseCallback("Bezier Curve", mouse_handler, nullptr);

    int key = -1;
    while (key != 27 && key != 's') 
    {
        for (auto &point : control_points) 
        {
            cv::circle(window, point, 1, {0,0,0}, 3);
        }

        if (control_points.size() == 4) 
        {
            naive_bezier(control_points, window);
            bezier(control_points, window);
            std::vector<cv::Point2f>().swap(control_points);    
        }

        cv::imshow("Bezier Curve", window);
        key = cv::waitKey(20);
    }
    if(key == 's')cv::imwrite("bezier_curve.ppm", window);
return 0;
}