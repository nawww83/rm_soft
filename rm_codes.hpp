#pragma once

#include <vector>
#include <span>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <concepts>
#include <complex>

#include "generals.hpp"

/**
 * @file rm_codec.hpp
 * @brief Модуль кодирования и быстрого мягкого декодирования кодов Рида-Маллера (RM).
 * 
 * Содержит функции для рекурсивного кодирования по схеме Плоткина, извлечения 
 * информационной последовательности, а также алгоритм быстрого мягкого декодирования 
 * на основе логарифмических отношений правдоподобия (LLR).
 */

/**
 * @brief Вычисление знаковой функции (Signum) для вещественных чисел.
 * 
 * @tparam T Вещественный тип данных (std::floating_point).
 * @param[in] val Входное значение.
 * @return T Возвращает 1.0, если val > 0; -1.0, если val < 0; иначе 0.0.
 */
template <std::floating_point T>
inline constexpr T sign(T val) {
    return (val > T(0)) - (val < T(0));
}

/**
 * @brief Внутренняя рекурсивная функция быстрого мягкого декодирования кода Рида-Маллера RM(r, m).
 * 
 * Алгоритм использует рекурсивное разбиение Плоткина \f$(u, u \oplus v)\f$ в терминах LLR.
 * Базовыми случаями рекурсии являются коды повторения RM(0, m) и коды с проверкой на четность RM(m-1, m).
 * Результат возвращается в биполярном формате (+1/-1).
 * 
 * @tparam T Вещественный тип данных для вычисления LLR.
 * @param[in] r Текущий порядок кода Рида-Маллера.
 * @param[in] m Текущий размер пространства (длина блока \f$N = 2^m\f$).
 * @param[in] llr Набор входных мягких решений (LLR) длины \f$2^m\f$.
 * @param[out] output Буфер для записи биполярного декодированного слова (+1/-1).
 * @param[in,out] ws Аллокатор памяти (Workspace) для временных векторов LLR на разных уровнях рекурсии.
 */
template <std::floating_point T>
inline void rm_soft_decode_recursive_fast(int r, int m, std::span<const T> llr, std::span<T> output, DecoderWorkspace<T>& ws) {
    const size_t n = llr.size();

    // Базовый случай 1: Код RM(0, m) — код с полным повторением.
    // Решение принимается по знаку суммы всех элементов LLR.
    if (r == 0) { 
        T sum = std::accumulate(llr.begin(), llr.end(), T(0));
        std::fill(output.begin(), output.end(), (sum >= T(0)) ? T(1) : T(-1));
        return;
    }
    
    // Базовый случай 2: Код RM(m-1, m) — код с одиночной проверкой на четность (SPC).
    // Используется мягкое декодирование по алгоритму Вагнера.
    if (r == m - 1) { 
        // Жесткое покомпонентное решение
        std::transform(llr.begin(), llr.end(), output.begin(), [](T x) { return (x >= T(0)) ? T(1) : T(-1); });
        // Проверка общего паритета произведения знаков
        T parity = std::reduce(output.begin(), output.end(), T(1), std::multiplies<T>());
        if (parity < T(0)) {
            // Если четность нарушена, инвертируем бит с наименьшей абсолютной величиной LLR
            auto min_it = std::min_element(llr.begin(), llr.end(), [](T a, T b) { return std::abs(a) < std::abs(b); });
            output[std::distance(llr.begin(), min_it)] *= T(-1);
        }
        return;
    }

    const size_t half = n / 2;
    auto llr1 = llr.subspan(0, half);
    auto llr2 = llr.subspan(half, half);

    // Шаг 1: Декодирование компоненты v2 (соответствует коду RM(r-1, m-1))
    auto llr_v2 = ws.allocate(half);
    auto v2_dec = ws.allocate(half);

    // Вычисление LLR для v2 по правилу f-функции (знак-минимум, аналогично полярным кодам)
    for (size_t i = 0; i < half; ++i) {
        llr_v2[i] = sign(llr1[i]) * sign(llr2[i]) * std::min(std::abs(llr1[i]), std::abs(llr2[i]));
    }
    rm_soft_decode_recursive_fast<T>(r - 1, m - 1, llr_v2, v2_dec, ws);

    // Шаг 2: Декодирование компоненты v1 (соответствует коду RM(r, m-1))
    auto llr_v1 = ws.allocate(half);
    auto v1_dec = ws.allocate(half);

    // Обновление LLR для v1 по правилу g-функции с использованием уже декодированного v2
    for (size_t i = 0; i < half; ++i) {
        llr_v1[i] = llr1[i] + llr2[i] * v2_dec[i];
    }
    rm_soft_decode_recursive_fast<T>(r, m - 1, llr_v1, v1_dec, ws);

    // Шаг 3: Объединение результатов по схеме Плоткина в биполярном домене
    for (size_t i = 0; i < half; ++i) {
        output[i] = v1_dec[i];
        output[i + half] = v1_dec[i] * v2_dec[i]; // Эквивалентно XOR в булевом домене
    }

    // Освобождение временной памяти в структуре Workspace
    ws.free(half * 4);
}

/**
 * @brief Интерфейсная функция быстрого мягкого декодирования кода Рида-Маллера RM(r, m).
 * 
 * Принимает вектор вещественных LLR, запускает рекурсивную процедуру и 
 * переводит итоговое кодовое слово из биполярного представления (+1/-1) в стандартные биты (0/1).
 * 
 * @tparam T Вещественный тип данных для LLR.
 * @param[in] r Порядок кода Рида-Маллера.
 * @param[in] m Параметр длины кода (\f$N = 2^m\f$).
 * @param[in] llr Мягкие решения из канала (логарифмические отношения правдоподобия).
 * @param[in,out] ws Аллокатор рабочей памяти для декодера.
 * @return std::vector<int> Декодированное кодовое слово в бинарном формате (0 или 1).
 */
template <std::floating_point T>
inline std::vector<int> rm_soft_decode_fast(int r, int m, std::span<const T> llr, DecoderWorkspace<T>& ws) {
    std::vector<T> bipolar_output(llr.size());
    ws.current_offset = 0;
    rm_soft_decode_recursive_fast<T>(r, m, llr, bipolar_output, ws);
    
    std::vector<int> decoded_bits(llr.size());
    for (size_t i = 0; i < llr.size(); ++i) {
        decoded_bits[i] = (bipolar_output[i] > 0) ? 0 : 1;
    }
    return decoded_bits;
}

/**
 * @brief Кодирование информационной последовательности кодом Рида-Маллера RM(r, m).
 * 
 * Функция вычисляет длину результирующего блока \f$N = 2^m\f$ и вызывает 
 * внутреннюю процедуру кодирования на базе рекурсивной конструкции Плоткина.
 * 
 * @param[in] r Порядок кода.
 * @param[in] m Параметр длины кодового слова.
 * @param[in] info Вектор информационных битов длины \f$K = \sum_{i=0}^{r} \binom{m}{i}\f$.
 * @return std::vector<int> Сгенерированное разрешенное кодовое слово длины \f$N\f$.
 */
inline std::vector<int> rm_encode(int r, int m, const std::vector<int>& info) {
    size_t n = 1 << m;
    std::vector<int> codeword(n);
    plotkin_encode_recursive(r, m, info, codeword);
    return codeword;
}

/**
 * @brief Рекурсивное извлечение исходных информационных битов из безошибочного кодового слова.
 * 
 * Процедура производит обратное разбиение Плоткина \f$(u, u \oplus v) \to (u, v)\f$ для 
 * систематического восстановления исходного информационного вектора.
 * 
 * @param[in] r Текущий порядок кода.
 * @param[in] m Текущий параметр длины блока.
 * @param[in] codeword Текущий сегмент кодового слова (или его промежуточная компонента).
 * @param[out] info Итоговый буфер, куда послойно записываются извлеченные информационные биты.
 */
inline void rm_extract_info_recursive(int r, int m, std::span<const int> codeword, std::span<int> info) {
    if (r == 0) { info[0] = codeword[0]; return; }
    if (r == m) { std::copy(codeword.begin(), codeword.end(), info.begin()); return; }

    const size_t half = codeword.size() / 2;
    int k_v1 = get_rm_k(r, m - 1);
    int k_v2 = get_rm_k(r - 1, m - 1);

    auto info_v1 = info.subspan(0, k_v1);
    auto info_v2 = info.subspan(k_v1, k_v2);

    std::vector<int> v1(codeword.begin(), codeword.begin() + half);
    std::vector<int> v2(half);
    // Обратное преобразование Плоткина в булевом домене: v = c1 ^ c2
    for (size_t i = 0; i < half; ++i) {
        v2[i] = codeword[i] ^ codeword[half + i];
    }

    rm_extract_info_recursive(r, m - 1, v1, info_v1);
    rm_extract_info_recursive(r - 1, m - 1, v2, info_v2);
}

/**
 * @brief Интерфейсная функция для извлечения информационных битов из кодового слова RM(r, m).
 * 
 * Выделяет память под вектор информационных битов требуемого размера \f$K\f$ 
 * и запускает рекурсивный алгоритм декомпозиции. Применяется после успешного жесткого декодирования.
 * 
 * @param[in] r Порядок кода Рида-Маллера.
 * @param[in] m Параметр длины кода.
 * @param[in] codeword Входное полностью декодированное ковое слово (длины \f$2^m\f$).
 * @return std::vector<int> Извлеченный вектор полезных информационных битов (длины \f$K\f$).
 */
inline std::vector<int> rm_extract_info(int r, int m, const std::vector<int>& codeword) {
    std::vector<int> info(get_rm_k(r, m));
    rm_extract_info_recursive(r, m, codeword, info);
    return info;
}
