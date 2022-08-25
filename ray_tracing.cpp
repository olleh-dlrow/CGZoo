/*
c++ ray_tracing.cpp -std=c++17 -o ray_tracing && ./ray_tracing
*/
#include<iostream>
#include<mutex>
using namespace std;

template<typename T>
struct Singleton {
    Singleton& getInstance() {
        static Singleton instance = Singleton();
        return instance;
    }

    // 目的：在ins创建完成之后，不要再用锁
    // 写法，空锁空
    Singleton& getIns() {
        {
            // can't put here, speed is low
            if(ins == nullptr) {
                unique_lock<mutex> lk(mt);
                if(ins == nullptr) {
                    ins = new Singleton();
                }
            }
        }

        return ins;
    }

    static Singleton* ins = nullptr;
    static mutex mt;
};

class A {
public:
};

class B: public A {
public:
    ~B() {}
    int b;
};

int main() {
    sizeof(B);
    cout << "ray tracing!\n";
    return 0;
}