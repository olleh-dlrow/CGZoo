/*
c++ hermite_curve.cpp -std=c++11 -o hermite_curve -lopencv_core -lopencv_highgui -lopencv_imgcodecs -lopencv_imgproc && ./hermite_curve
*/
#include<opencv2/opencv.hpp>
#include<iostream>
using std::vector;

vector<cv::Point2d> control_points;

void mouse_handler(int event, int x, int y, int flags, void *userdata) 
{
    if (event == cv::EVENT_LBUTTONDOWN) 
    {
        std::cout << "Left button of the mouse is clicked - position (" << x << ", "
        << y << ")" << '\n';
        control_points.emplace_back(x, y);
    }     
}

void hermite_curve(const vector<cv::Point2d>& control_points, cv::Mat& window) {
    cv::Point2d P = control_points[1];
    cv::Point2d Q = control_points[2];
    cv::Vec2d   v = control_points[1] - control_points[0];
    cv::Vec2d   w = control_points[3] - control_points[2];
    cv::Mat G = (cv::Mat_<double>(2, 4) << 
                P.x, Q.x, v[0], w[0],
                P.y, Q.y, v[1], w[1]);
    cv::Mat M = (cv::Mat_<double>(4, 4) << 
                1, 0, -3, 2, 
                0, 0,  3,-2,
                0, 1, -2, 1,
                0, 0, -1, 1);
    double eps = 1e-3;
    for(double t = 0.0; t <= 1.0; t += eps) {
        cv::Mat T = (cv::Mat_<double>(4, 1) << 1, t, t*t, t*t*t);
        cv::Mat gamma = G * M * T;  // 2 * 1
        double x = gamma.at<double>(0, 0);
        double y = gamma.at<double>(1, 0);
        window.at<cv::Vec3b>(int(y), int(x)) = cv::Vec3b(255, 0, 0);           
    }       
}

int main() {
    cv::Mat window = cv::Mat(700, 700, CV_8UC3, cv::Vec3b(255, 255, 255));
    cv::cvtColor(window, window, cv::COLOR_BGR2RGB);        // the method of storing is BGR
    cv::namedWindow("Hermite Curve", cv::WINDOW_AUTOSIZE);  // you can't change the size of window
    cv::setMouseCallback("Hermite Curve", mouse_handler, nullptr); 

    int key = -1;
    while (key != 27 && key != 's')
    {
        for (auto &point : control_points) 
        {
            cv::circle(window, point, 1, {0,0,0}, 3);
        }

        if (control_points.size() == 4) 
        {
            hermite_curve(control_points, window);
            std::vector<cv::Point2d>().swap(control_points);    
        }

        cv::imshow("Hermite Curve", window);
        key = cv::waitKey(20);
    }

    if(key == 's')cv::imwrite("hermite_curve.ppm", window);
    return 0;
}