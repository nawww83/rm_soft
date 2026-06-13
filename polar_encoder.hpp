#include <vector>
#include <span>


/**
 * @file polar_encoder.hpp
 * @brief Блок рекурсивного полярного кодирования (GF(2) преобразование).
 */

/**
 * @brief Рекурсивное инкрементальное ядро кодирования. Реализует операцию "бабочки".
 * @param[in,out] codeword Спан данных, трансформируемый на месте из вектора u в кодовое слово x.
 */
inline void polar_encode_recursive_core(std::span<int> codeword) {
    const size_t n = codeword.size();
    if (n <= 1) return;
    const size_t half = n / 2;
    
    // x1 = u1 ^ u2, x2 = u2
    for (size_t i = 0; i < half; ++i) {
        codeword[i] = codeword[i] ^ codeword[half + i];
    }
    
    // Рекурсивный спуск к нижним подблокам
    polar_encode_recursive_core(codeword.subspan(0, half));
    polar_encode_recursive_core(codeword.subspan(half, half));
}

/**
 * @brief Интерфейсная функция кодирования блока данных.
 * 
 * @param n Целевая длина полярного блока.
 * @param info Чистый вектор информационных битов (длина K).
 * @param bit_mask Сгенерированная маска заморозки 5G.
 * @return std::vector<int> Выходное кодовое слово (длина N).
 */
inline std::vector<int> pure_polar_encode(size_t n, const std::vector<int>& info, const std::vector<uint8_t>& bit_mask) {
    std::vector<int> u(n, 0);
    size_t info_idx = 0;
    for (size_t i = 0; i < n; ++i) {
        if (bit_mask[i] == 1) {
            u[i] = info[info_idx++];
        }
    }
    std::vector<int> codeword = u;
    polar_encode_recursive_core(codeword);
    return codeword;
}
