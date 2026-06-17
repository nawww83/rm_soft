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

inline void plotkin_encode_recursive(int r, int m, std::span<const int> info, std::span<int> output, DecoderWorkspace<int>& ws) {
    const size_t n = output.size();
    
    if (r == 0) { 
        assert(!info.empty() && "Info span cannot be empty in base case r=0");
        std::fill(output.begin(), output.end(), info[0]); 
        return; 
    }
    if (r == m) { 
        assert(info.size() == n && "Info size must match output size in base case r=m");
        std::copy(info.begin(), info.end(), output.begin()); 
        return; 
    }

    int k_v1 = get_rm_k(r, m - 1);
    int k_v2 = get_rm_k(r - 1, m - 1);

    assert(info.size() == static_cast<size_t>(k_v1 + k_v2) && "Total info size mismatch with k_v1 + k_v2");

    auto info_v1 = info.subspan(0, k_v1);
    auto info_v2 = info.subspan(k_v1, k_v2);

    const size_t half = n / 2;

    typename DecoderWorkspace<int>::Guard guard(ws);
    std::span<int> v1 = ws.allocate(half);
    std::span<int> v2 = ws.allocate(half);

    plotkin_encode_recursive(r, m - 1, info_v1, v1, ws);
    plotkin_encode_recursive(r - 1, m - 1, info_v2, v2, ws);

    for (size_t i = 0; i < half; ++i) {
        output[i] = v1[i];
        output[i + half] = v1[i] ^ v2[i];
    }
}
