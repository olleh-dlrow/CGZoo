// https://www.zhihu.com/question/27908489/answer/355105668

/*
c++ multithread_test.cpp -std=c++17 -o multithread_test && ./multithread_test
*/

#include <mutex>
#include <condition_variable>
#include <functional>
#include <queue>
#include <thread>
#include <iostream>
#include <unistd.h>

class fixed_thread_pool {
 public:
  explicit fixed_thread_pool(size_t thread_count)
      : data_(std::make_shared<data>()) {
    for (size_t i = 0; i < thread_count; ++i) {
      ths.emplace_back([data = data_] {
        std::unique_lock<std::mutex> lk(data->mtx_);
        for (;;) {
          if (!data->tasks_.empty()) {
            auto current = std::move(data->tasks_.front());
            data->tasks_.pop();
            lk.unlock();
            current();
            lk.lock();
          } else if (data->is_shutdown_) {
            break;
          } else {
            data->cond_.wait(lk);
          }
        }
      });
    }
  }

  fixed_thread_pool() = default;
  fixed_thread_pool(fixed_thread_pool&&) = default;

  ~fixed_thread_pool() {
    if ((bool) data_) {
      {
        std::lock_guard<std::mutex> lk(data_->mtx_);
        data_->is_shutdown_ = true;
      }
      data_->cond_.notify_all();
      for (auto &&t : ths)
      {
        t.join();
      }
      
    }
  }

  template <class F>
  void execute(F&& task) {
    {
      std::lock_guard<std::mutex> lk(data_->mtx_);
      data_->tasks_.emplace(std::forward<F>(task));
    }
    data_->cond_.notify_one();
  }

 private:
 std::vector<std::thread> ths;
  struct data {
    std::mutex mtx_;
    std::condition_variable cond_;
    bool is_shutdown_ = false;
    std::queue<std::function<void()>> tasks_;
  };
  std::shared_ptr<data> data_;
};

const int W = 10, H = 10;

int pixels[W*H];
int mm = 0;
int main() {
  {
    fixed_thread_pool pool(10);
    for(int i = 0; i < W*H; i++) {
        pool.execute([i](){
            pixels[i] = i;
            std::cout << "sleep " << i << std::endl;
            sleep(1);
        });
    }
  }
    std::copy(pixels, pixels + W*H, std::ostream_iterator<int>(std::cout, " "));
    return 0;
}