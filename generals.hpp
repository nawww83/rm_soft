#pragma once

#include <vector>
#include <span>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <concepts>
#include <complex>
#include <cassert>

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
        // Выделяем память с хорошим запасом под дерево LLR и решений
        memory_pool.resize(max_size * 8); 
    }

    std::span<T> allocate(size_t size) {
        // Защита от выхода за пределы пула
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
            // При уничтожении гарда offset гарантированно возвращается в исходное состояние
            ws.current_offset = saved_offset;
        }
    };
};
