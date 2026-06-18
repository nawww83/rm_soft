#pragma once

#include <vector>
#include <span>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <concepts>
#include <complex>
#include <cassert>
#include <future>
#include <thread>
#include <queue>
#include <mutex>
#include <condition_variable>
#include <functional>

#include "nr_5g_polar_table.hpp"

enum class CodeType
{
    ReedMuller,
    PolarFastSsc,
    PolarSc
};

class UniversalRandomInterleaver
{
private:
    size_t N;
    std::vector<size_t> p_forward;
    std::vector<size_t> p_backward;

public:
    explicit UniversalRandomInterleaver(size_t block_size, unsigned int seed = 42) : N(block_size)
    {
        p_forward.resize(N);
        std::iota(p_forward.begin(), p_forward.end(), 0);
        std::mt19937 g(seed);
        std::shuffle(p_forward.begin(), p_forward.end(), g);

        p_backward.resize(N);
        for (size_t i = 0; i < N; ++i)
        {
            p_backward[p_forward[i]] = i;
        }
    }

    size_t size() const { return N; }

    inline void interleave(std::span<const int> src, std::span<int> dst) const
    {
        for (size_t i = 0; i < N; ++i)
        {
            dst[i] = src[p_forward[i]];
        }
    }

    template <std::floating_point T>
    inline void deinterleave(std::span<const T> src, std::span<T> dst) const
    {
        for (size_t i = 0; i < N; ++i)
        {
            dst[i] = src[p_backward[i]];
        }
    }
};

inline std::vector<uint8_t> get_polar_5g_mask(size_t n, int k) {
    assert(n <= 1024 && (n & (n - 1)) == 0 && "N должно быть степенью 2");
    std::vector<int> filtered;
    filtered.reserve(n);
    for (uint16_t idx : nr_5g_polar::reliability_table_1024) {
        if (idx < static_cast<int>(n)) {
            filtered.push_back(idx);
        }
    }
    std::vector<uint8_t> channel_mask(n, 0);
    size_t start_info_idx = filtered.size() - k;
    for (size_t i = start_info_idx; i < filtered.size(); ++i) {
        channel_mask[filtered[i]] = 1;
    }
    return channel_mask;
}

inline constexpr int binomial(int n, int k) {
    if (k < 0 || k > n) return 0;
    if (k == 0 || k == n) return 1;
    if (k > n / 2) k = n - k;
    int res = 1;
    for (int i = 1; i <= k; ++i) res = res * (n - i + 1) / i;
    return res;
}

inline constexpr int get_rm_k(int r, int m) {
    int k = 0;
    for (int i = 0; i <= r; ++i) k += binomial(m, i);
    return k;
}

template <typename T>
struct DecoderWorkspace {
    std::vector<T> memory_pool;
    size_t current_offset = 0;

    explicit DecoderWorkspace(size_t max_size) {
        memory_pool.resize(max_size * 8); 
    }

    std::span<T> allocate(size_t size) {
        if (current_offset + size > memory_pool.size()) {
            std::cerr << "CRITICAL: Workspace pool overflow! Offset: " << current_offset 
                      << ", Requested: " << size << ", Capacity: " << memory_pool.size() << std::endl;
            assert(false);
        }
        std::span<T> sp(memory_pool.data() + current_offset, size);
        current_offset += size;
        return sp;
    }

    // Класс для автоматического отката offset (RAII Стек)
    struct Guard {
        DecoderWorkspace& ws;
        size_t saved_offset;

        explicit Guard(DecoderWorkspace& workspace) 
            : ws(workspace), saved_offset(workspace.current_offset) {}

        ~Guard() {
            ws.current_offset = saved_offset;
        }
    };
};

// Простая реализация пула потоков (Thread Pool)
class ThreadPool
{
private:
    std::vector<std::jthread> workers;
    std::queue<std::function<void()>> tasks;
    std::mutex queue_mutex;
    std::condition_variable cv;
    bool stop = false;

public:
    ThreadPool(size_t threads)
    {
        for (size_t i = 0; i < threads; ++i)
        {
            workers.emplace_back([this](std::stop_token st)
                                 {
                while (!st.stop_requested()) {
                    std::function<void()> task;
                    {
                        std::unique_lock<std::mutex> lock(this->queue_mutex);
                        this->cv.wait(lock, [this, &st] { 
                            return this->stop || !this->tasks.empty() || st.stop_requested(); 
                        });
                        
                        if ((this->stop || st.stop_requested()) && this->tasks.empty())
                            return;
                            
                        task = std::move(this->tasks.front());
                        this->tasks.pop();
                    }
                    task(); // Выполняем симуляцию точки
                } });
        }
    }

    template <class F, class... Args>
    auto enqueue(F &&f, Args &&...args)
        -> std::future<typename std::invoke_result<F, Args...>::type>
    {
        using return_type = typename std::invoke_result<F, Args...>::type;

        auto task = std::make_shared<std::packaged_task<return_type()>>(
            std::bind(std::forward<F>(f), std::forward<Args>(args)...));

        std::future<return_type> res = task->get_future();
        {
            std::lock_guard<std::mutex> lock(queue_mutex);
            if (stop)
                throw std::runtime_error("Добавление задачи в остановленный пул");
            tasks.emplace([task]()
                          { (*task)(); });
        }
        cv.notify_one();
        return res;
    }

    ~ThreadPool()
    {
        {
            std::lock_guard<std::mutex> lock(queue_mutex);
            stop = true;
        }
        cv.notify_all();
        // jthread автоматически сделают join при уничтожении вектора workers
    }
};