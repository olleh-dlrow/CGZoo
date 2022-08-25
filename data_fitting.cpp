/*
c++ data_fitting.cpp -std=c++11 -o data_fitting -lopencv_core -lopencv_highgui -lopencv_imgcodecs -lopencv_imgproc && ./data_fitting
*/
#include <iostream>
#include <opencv2/opencv.hpp>

#define CVUI_IMPLEMENTATION
#include "cvui.h"

const char* VIEWPORT_NAME = "Viewport";
const char* CONSOLE_NAME  = "Console";
const char* SAVE_FILE     = "data_fitting.ppm";

const cv::Vec3b RED     = cv::Vec3b(0, 0, 255);
const cv::Vec3b GREEN   = cv::Vec3b(0, 255, 0);
const cv::Vec3b BLUE    = cv::Vec3b(255, 0, 0);
const cv::Vec3b PURPLE  = cv::Vec3b(255, 0, 255);
const cv::Vec3b BLACK   = cv::Vec3b(0, 0, 0);
const cv::Vec3b WHITE   = cv::Vec3b(255, 255, 255);
const cv::Vec3b BG      = BLACK;

const int CTL_POINT = 0;
const int CTL_LINE  = 1;
int ctl_mode = CTL_POINT;

std::vector<cv::Point2f> control_points;
cv::Mat window;
cv::Mat console_window;
int m = 5;
float d_m = 1;

float lambda = 0.01;
float d_lambda = 0.01;

float sigma = 1.0;
float d_sigma = 0.1;

// float* var  [] = {&m, &lambda, &sigma};
// float* d_var[] = {&d_m, &d_lambda, &d_sigma};
// int    var_i = 0;

void draw(cv::Mat& window, int x, int y, cv::Vec3b color) { if(x >= 0 && x < window.cols && y >= 0 && y < window.cols)window.at<cv::Vec3b>(y, x) = color;   }
void mouse_handler(int event, int x, int y, int flags, void *userdata) 
{
    if (event == cv::EVENT_LBUTTONDOWN) 
    {
        // std::cout << "Left button of the mouse is clicked - position (" << x << ", " << y << ")" << '\n';
        control_points.emplace_back(x, y);
    }     
}
void polynomial_interpolation_fitting(const std::vector<cv::Point2f>& points, cv::Mat& window) {
    int     n     = points.size() - 1;
    cv::Mat b     = cv::Mat(n + 1, 1, CV_32F);
    cv::Mat A     = cv::Mat(n + 1, n + 1, CV_32F);
    cv::Mat alpha = cv::Mat(n + 1, 1, CV_32F);

    for(int i = 0; i <= n; i++) {
        b.at<float>(i) = points[i].y;
    }
    for(int i = 0; i <= n; i++) {
        for(int j = 0; j <= n; j++) {
            A.at<float>(i, j) = std::pow(points[i].x, j);
        }
    }
    alpha = A.inv() * b;

    float x0  = points[0].x;
    float xn  = points[n].x;
    float eps = 1e-1;

    for(float x = x0; x <= xn; x += eps) {
        float y = 0;
        for(int i = 0; i <= n; i++) {
            y += alpha.at<float>(i) * std::pow(x, i);
        }
        draw(window, x, y, BLUE);
    }
}
void Gauss_interpolation_fitting(const std::vector<cv::Point2f>& points, cv::Mat& window) {
    int     n     = points.size() - 1;
    cv::Mat b     = cv::Mat(n + 1, 1, CV_32F);
    cv::Mat A     = cv::Mat(n + 1, n + 1, CV_32F);
    cv::Mat alpha = cv::Mat(n + 1, 1, CV_32F);
    float b0  = 0;

    for(int i = 0; i <= n; i++) {
        b.at<float>(i) = points[i].y;
    }
    b0 = cv::mean(b)[0];
    for(int i = 0; i <= n; i++) {
        b.at<float>(i) -= b0;
    }

    for(int i = 0; i <= n; i++) {
        for(int j = 0; j <= n; j++) {
            A.at<float>(i, j) = std::exp(-1.0 / (2.0 * sigma*sigma) * std::pow((points[i].x - points[j].x), 2));
        }
    }
    alpha = A.inv() * b;  

    float x0  = points[0].x;
    float xn  = points[n].x;
    float eps = 1e-1;

    for(float x = x0; x <= xn; x += eps) {
        float y = b0;
        for(int i = 0; i <= n; i++) {
            y += alpha.at<float>(i) * std::exp(-1.0 /  (2.0 * sigma*sigma) * std::pow((x - points[i].x), 2));
        }
        draw(window, x, y, GREEN);
    }          
}
void polynomial_regression(const std::vector<cv::Point2f>& points, cv::Mat& window, int m) {
    int n = points.size() - 1;
    if(m >= n + 1)return;
    cv::Mat X = cv::Mat(n + 1, m + 1, CV_32F);
    cv::Mat Y = cv::Mat(n + 1, 1, CV_32F);
    cv::Mat W = cv::Mat(m + 1, 1, CV_32F);

    for(int i = 0; i <= n; i++) {
        for(int j = 0; j <= m; j++) {
            X.at<float>(i, j) = std::pow(points[i].x, j);
        }
    }
    for(int i = 0; i <= n; i++) {
        Y.at<float>(i) = points[i].y;
    }
    W = (X.t() * X).inv() * X.t() * Y;

    float x0  = points[0].x;
    float xn  = points[n].x;
    float eps = 1e-1;

    for(float x = x0; x <= xn; x += eps) {
        float y = 0;
        for(int i = 0; i <= m; i++) {
            y += W.at<float>(i) * std::pow(x, i);
        }
        draw(window, x, y, RED);
    }
}
void ridge_regression(const std::vector<cv::Point2f>& points, cv::Mat& window, int m, float lambda) {
    int n = points.size() - 1;
    if(m >= n + 1)return;
    cv::Mat X = cv::Mat(n + 1, m + 1, CV_32F);
    cv::Mat Y = cv::Mat(n + 1, 1, CV_32F);
    cv::Mat W = cv::Mat(m + 1, 1, CV_32F);
    cv::Mat I = cv::Mat::eye(m + 1, m + 1, CV_32F);

    for(int i = 0; i <= n; i++) {
        for(int j = 0; j <= m; j++) {
            X.at<float>(i, j) = std::pow(points[i].x, j);
        }
    }
    for(int i = 0; i <= n; i++) {
        Y.at<float>(i) = points[i].y;
    }
    W = (X.t() * X + lambda * I).inv() * X.t() * Y;

    float x0  = points[0].x;
    float xn  = points[n].x;
    float eps = 1e-1;

    for(float x = x0; x <= xn; x += eps) {
        float y = 0;
        for(int i = 0; i <= m; i++) {
            y += W.at<float>(i) * std::pow(x, i);
        }
        draw(window, x, y, PURPLE);
    }
}
void plot() {
    polynomial_interpolation_fitting(control_points, window);
    Gauss_interpolation_fitting(control_points, window);
    polynomial_regression(control_points, window, m);
    ridge_regression(control_points, window, m, lambda);
}
void clear() {  window = cv::Mat(700, 700, CV_8UC3, BG);    }
void print_var() {
    printf("m:%d  lambda:%.3f  sigma:%.3f\t\n", (int)m, lambda, sigma); 
}
void replot() {
    // print_var();
    clear();
    plot();
}

void handle_viewport() {
    cvui::context(VIEWPORT_NAME);

    if(ctl_mode == CTL_POINT) {
        if(cvui::mouse(cvui::CLICK)) {
            int x = cvui::mouse().x;
            int y = cvui::mouse().y;
            control_points.emplace_back(x, y);
            printf("clicked\n");
        }
        for (auto &point : control_points) {
            cv::circle(window, point, 1, WHITE, 3);
        }
    }
    if(ctl_mode == CTL_LINE) {
        replot();
    }

    cvui::imshow(VIEWPORT_NAME, window);
}

void handle_console() {
    console_window = cv::Mat(700, 300, CV_8UC3, cv::Scalar(49, 52, 49));
    cvui::context(CONSOLE_NAME);

    // cvui::trackbar(250, &sigma, 0.1f, 20.0f, 1, "%.2Lf");
    cvui::beginColumn(console_window, 10, 50, 300, 700, 10);
    cvui::trackbar(200, &sigma, 0.1f, 20.0f, 1, "%.2Lf");         // Gauss...
    cvui::trackbar(200, &m, 0, 10, 1, "%.0Lf");                    // regression polynomial...
    cvui::trackbar(200, &lambda, 0.0f, 10.0f, 1, "%.2Lf");       // regression punish
    if(cvui::button("CLEAR")) {
        clear();
        std::vector<cv::Point2f>().swap(control_points);
        ctl_mode = CTL_POINT;
    }
    if(ctl_mode == CTL_POINT && cvui::button("DRAW")) {
        ctl_mode = CTL_LINE;
    }
    cvui::endColumn();           

    cvui::imshow(CONSOLE_NAME, console_window);
}

int main() 
{
    int delay_time = 20;
    window = cv::Mat(700, 700, CV_8UC3, BG);
    
    const cv::String windows[] = {VIEWPORT_NAME, CONSOLE_NAME};
    cvui::init(windows, 2);

    while(true) {
        handle_console();
        handle_viewport();
        if(cv::waitKey(delay_time) == 27) {
            break;
        }
    }

    // cv::namedWindow(w_name, cv::WINDOW_AUTOSIZE);
    // cv::setMouseCallback(w_name, mouse_handler, nullptr);

    // int key = -1;
    // while (key != 27) 
    // {
    //     for (auto &point : control_points) {
    //         cv::circle(window, point, 1, WHITE, 3);
    //     }
    //     if(key >= '0' && key <= sizeof(var)/sizeof(float*) + '0' - 1) {
    //         var_i = key - '0';
    //         std::cout << "\nchanged to ctl " << var_i << "\n";
    //     }
    //     switch(key) {
    //     case 'h':
    //         std::cout << "\nhelp:\n0:m\t(+-1)\n1:lambda\t(+-0.01)\n2:sigma\t(+-0.1)\n";
    //         std::cout << "now at ctl " << var_i << "\n";
    //         break;
    //     case 't':
    //         std::cout << "\nstatus:\n";
    //         print_var();
    //         std::cout << "now at ctl " << var_i << "\n";
    //         break;
    //     case 'd':
    //         plot();
    //         std::cout << "\ndrawed\n";
    //         break;
    //     case 'c':
    //         clear();
    //         std::vector<cv::Point2f>().swap(control_points);
    //         std::cout << "\ncleared\n";
    //         break;
    //     case 's':
    //         cv::imwrite(save_file, window);
    //         std::cout << "\nsaved\n";
    //         break;
    //     case '[':
    //         *var[var_i] -= *d_var[var_i];
    //         replot();
    //         break;
    //     case ']':
    //         *var[var_i] += *d_var[var_i];
    //         replot();
    //         break;
    //     default:break;
    //     }
    //     cv::imshow(w_name, window);
    //     key = cv::waitKey(20);
    // }
    return 0;
}
