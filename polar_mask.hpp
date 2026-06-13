#pragma once

#include <vector>
#include <cstdint>
#include <cmath>
#include <numeric>
#include <algorithm>

namespace polar_mask {

// Аппроксимация функции Гаусса для полярных кодов
inline double phi(double x) {
    if (x < 0) return 1.0;
    if (x <= 10.0) {
        return std::exp(-0.4527 * std::pow(x, 0.86) + 0.0218);
    } else {
        return std::sqrt(M_PI / x) * std::exp(-x / 4.0) * (1.0 - 10.0 / (7.0 * x));
    }
}

// Обратная функция phi^(-1)
inline double phi_inv(double y) {
    if (y >= 0.999) return 0.0001;
    if (y <= 0.0001) return 50.0;
    
    double low = 0.0001, high = 50.0;
    for (int iter = 0; iter < 20; ++iter) {
        double mid = low + (high - low) / 2.0;
        if (phi(mid) < y) high = mid; // phi монотонно убывает
        else low = mid;
    }
    return low;
}


inline std::vector<uint8_t> generate_custom_qam_mask(size_t n, int k, double eb_n0_db) {
    // 1. Рассчитываем базовый уровень энергии для текущего Eb/N0
    double snr = std::pow(10.0, eb_n0_db / 10.0);
    
    // Средние LLR для «сильных» и «слабых» бит 16-QAM
    double m_high = 2.0 * snr * (4.0 / std::sqrt(10.0)); // Коэффициенты нормировки созвездия 3GPP
    double m_low  = 2.0 * snr * (2.0 / std::sqrt(10.0));
    
     // ВАЖНО: Интерливер усреднил канал! 
    // Все N входных каналов получают одинаковую стартовую надежность.
    double m_avg = (m_high + m_low) / 2.0;
    
    std::vector<double> m(n, m_avg); // Заполняем весь вектор одинаковыми m_avg

    // 2. Эволюция плотности сверху вниз по уровням дерева (Кронекеровский обход)
    size_t steps = std::bit_width(n) - 1; // log2(N)
    for (size_t step = 0; step < steps; ++step) {
        size_t block_size = 1 << (step + 1);
        size_t half = block_size / 2;
        
        for (size_t i = 0; i < n; i += block_size) {
            for (size_t j = 0; j < half; ++j) {
                double m1 = m[i + j];
                double m2 = m[i + half + j];
                
                // Эволюция для левой ветви (операция f / бабочка)
                double m_left = phi_inv(1.0 - (1.0 - phi(m1)) * (1.0 - phi(m2)));
                
                // Эволюция для правой ветви (операция g)
                double m_right = m1 + m2;
                
                m[i + j] = m_left;
                m[i + half + j] = m_right;
            }
        }
    }

    // 3. Сортировка полученных каналов по надежности
    // Создаем пары: {индекс_канала, его_финальная_надежность}
    std::vector<std::pair<size_t, double>> reliability;
    for (size_t i = 0; i < n; ++i) {
        reliability.push_back({i, m[i]});
    }
    
    // Сортируем так, чтобы самые надежные (с максимальным m) оказались в конце
    std::sort(reliability.begin(), reliability.end(), [](const auto& a, const auto& b) {
        return a.second < b.second;
    });

    // 4. Формируем маску: K лучшим позициям даем 1, остальным 0
    std::vector<uint8_t> mask(n, 0);
    for (size_t i = n - k; i < n; ++i) {
        mask[reliability[i].first] = 1;
    }

    return mask;
}

}