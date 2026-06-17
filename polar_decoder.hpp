#pragma once

#include <vector>
#include <span>
#include <cmath>
#include <algorithm>
#include <concepts>
#include "generals.hpp"
#include "node_types.hpp"

/**
 * @file polar_decoder.hpp
 * @brief Блок быстрого мягкого Fast-SSC декодирования полярных кодов стандарта 5G NR.
 *
 * Логика вычислений f и g функций строго соответствует спецификации 3GPP TS 38.212,
 * где положительный LLR соответствует логическому 0, а отрицательный — логической 1.
 */

/**
 * @brief Внутренняя рекурсивная функция декодирования по полярному дереву (Fast-SSC).
 *
 * Выполняет последовательное исключение (Successive Cancellation) с оптимизацией
 * узлов типа Rate-0, Rate-1 и Repetition (REP) без спуска к листьям.
 *
 * @tparam T Вещественный тип вычислений LLR (float / double).
 * @param leaf_offset Индекс смещения текущего узла относительно начала глобальной маски.
 * @param llr Входной массив LLR для текущего подблока.
 * @param output_codeword Выходной массив для записи биполярного кодового подслова (+1/-1).
 * @param bit_mask Ссылка на полярную маску 5G.
 * @param ws Ссылка на контекст арены памяти (Workspace).
 */
template <std::floating_point T>
inline void fast_ssc_recursive(size_t leaf_offset, std::span<const T> llr,
                                        std::span<T> output_codeword,
                                        const std::vector<uint8_t> &bit_mask, DecoderWorkspace<T> &ws)
{
    const size_t n = llr.size();
    NodeType type = determine_node_type(leaf_offset, n, bit_mask);

    // =========================================================================
    // БЫСТРЫЕ УЗЛЫ ОСТАНОВА (Fast-SSC Оптимизации)
    // =========================================================================
    if (type == NodeType::Rate0)
    {
        // Все биты поддерева заморожены (u = 0). Согласно 3GPP, логический 0 отображается в +1.0
        std::fill(output_codeword.begin(), output_codeword.end(), T(1));
        return;
    }
    if (type == NodeType::Rate1)
    {
        // Узел чистой информации. Выполняется жесткое решение по знаку LLR
        for (size_t i = 0; i < n; ++i)
        {
            output_codeword[i] = (llr[i] >= T(0)) ? T(1) : T(-1);
        }
        return;
    }
    if (type == NodeType::Rep)
    {
        // Узел повторения. Решение принимается по знаку суммы всех LLR блока
        T sum_llr = 0;
        for (size_t i = 0; i < n; ++i)
            sum_llr += llr[i];
        T bit_sign = (sum_llr >= T(0)) ? T(1) : T(-1);
        std::fill(output_codeword.begin(), output_codeword.end(), bit_sign);
        return;
    }
    if (n == 1)
    {
        output_codeword[0] = (llr[0] >= T(0)) ? T(1) : T(-1);
        return;
    }

    // =========================================================================
    // КЛАССИЧЕСКИЙ СТАНДАРТНЫЙ СПУСК SC (Generic узел)
    // =========================================================================
    const size_t half = n / 2;
    typename DecoderWorkspace<T>::Guard memory_guard(ws); // Авто-откат смещения арены (RAII)

    auto llr_v1 = ws.allocate(half);
    auto v1_dec = ws.allocate(half);
    auto llr1 = llr.subspan(0, half);
    auto llr2 = llr.subspan(half, half);

    // Левая ветвь: f-операция (MIN-SUM аппроксимация "знак-минимум")
    for (size_t i = 0; i < half; ++i)
    {
        T s1 = (llr1[i] >= T(0)) ? T(1) : T(-1);
        T s2 = (llr2[i] >= T(0)) ? T(1) : T(-1);
        llr_v1[i] = s1 * s2 * std::min(std::abs(llr1[i]), std::abs(llr2[i]));
    }
    fast_ssc_recursive<T>(leaf_offset, llr_v1, v1_dec, bit_mask, ws);

    auto llr_v2 = ws.allocate(half);
    auto v2_dec = ws.allocate(half);

    // Правая ветвь: g-операция 3GPP (знак ПЛЮС согласован со спецификацией TS 38.212)
    for (size_t i = 0; i < half; ++i)
    {
        llr_v2[i] = llr2[i] + llr1[i] * v1_dec[i];
    }
    fast_ssc_recursive<T>(leaf_offset + half, llr_v2, v2_dec, bit_mask, ws);

    // Объединение решений (схлопывание узла Плоткина: x1 = v1 * v2, x2 = v2)
    for (size_t i = 0; i < half; ++i)
    {
        output_codeword[i] = v1_dec[i] * v2_dec[i];
        output_codeword[i + half] = v2_dec[i];
    }
}

/**
 * @brief Верхнеуровневый интерфейс декодера с прямым извлечением вектора u.
 *
 * @tparam T Вещественный тип вычислений LLR.
 * @param llr Входной массив LLR из канала.
 * @param bit_mask Ссылка на полярную маску 5G.
 * @param ws Контекст предвыделенной рабочей памяти.
 * @return std::vector<int> Декодированный вектор u (содержит инфо и замороженные биты).
 */
template <std::floating_point T>
inline std::vector<int> polar_decode_fast_generic(std::span<const T> llr, const std::vector<uint8_t> &bit_mask, DecoderWorkspace<T> &ws)
{
    const size_t n = llr.size();
    std::vector<T> bipolar_output(n);

    typename DecoderWorkspace<T>::Guard top_guard(ws);
    ws.current_offset = 0;

    // 1. Запускаем Fast-SSC рекурсию. bipolar_output теперь содержит очищенный вектор x.
    fast_ssc_recursive<T>(0, llr, bipolar_output, bit_mask, ws);

    // 2. Итеративный биполярный разворот (In-place декранчирование (разворот) матрицы Арикана)
    // Идем от мелких блоков к крупным, выполняя бабочки Плоткина наоборот.
    for (size_t len = 2; len <= n; len <<= 1)
    {
        size_t half = len / 2;
        for (size_t i = 0; i < n; i += len)
        {
            for (size_t j = 0; j < half; ++j)
            {
                T x1 = bipolar_output[i + j];
                T x2 = bipolar_output[i + half + j];

                // В биполярном представлении: u1 = x1 * x2, u2 = x2
                bipolar_output[i + j] = x1 * x2;
                // x2 остается на своем месте (bipolar_output[i + half + j] = x2), его не меняем
            }
        }
    }

    // 3. Перевод очищенного in-place вектора u в логические биты 0/1 для 5G
    std::vector<int> u(n);
    for (size_t i = 0; i < n; ++i)
    {
        u[i] = (bipolar_output[i] > 0) ? 0 : 1;
    }

    return u;
}

/**
 * @brief Честная рекурсивная функция SC-декодера.
 * Записывает декодированные информационные и замороженные биты прямо в вектор u.
 */
template <std::floating_point T>
inline void sc_recursive(size_t leaf_offset, std::span<const T> llr,
                                      std::span<T> output_codeword, // Нужен для g-операции соседа
                                      std::vector<int> &u_bits,       // Сюда пишем финальный ответ
                                      const std::vector<uint8_t> &bit_mask, DecoderWorkspace<T> &ws)
{
    const size_t n = llr.size();

    // Базовый случай: мы в листе дерева. Декодируем конкретный бит u[leaf_offset]
    if (n == 1)
    {
        if (bit_mask[leaf_offset] == 0) // Замороженный бит
        {
            u_bits[leaf_offset] = 0;    // Логический 0
            output_codeword[0] = T(1);  // Биполярный +1.0 для g-операции
        }
        else // Информационный бит
        {
            // Жесткое решение по знаку LLR: >= 0 это 0, < 0 это 1
            u_bits[leaf_offset] = (llr[0] >= T(0)) ? 0 : 1;
            output_codeword[0] = (llr[0] >= T(0)) ? T(1) : T(-1);
        }
        return;
    }

    const size_t half = n / 2;
    typename DecoderWorkspace<T>::Guard memory_guard(ws);

    auto llr_v1 = ws.allocate(half);
    auto v1_dec = ws.allocate(half);
    auto llr1 = llr.subspan(0, half);
    auto llr2 = llr.subspan(half, half);

    // Левая ветвь: f-операция
    for (size_t i = 0; i < half; ++i)
    {
        T s1 = (llr1[i] >= T(0)) ? T(1) : T(-1);
        T s2 = (llr2[i] >= T(0)) ? T(1) : T(-1);
        llr_v1[i] = s1 * s2 * std::min(std::abs(llr1[i]), std::abs(llr2[i]));
    }
    // Спускаемся влево
    sc_recursive<T>(leaf_offset, llr_v1, v1_dec, u_bits, bit_mask, ws);

    auto llr_v2 = ws.allocate(half);
    auto v2_dec = ws.allocate(half);

    // Правая ветвь: g-операция (использует v1_dec, полученный из левой ветви)
    for (size_t i = 0; i < half; ++i)
    {
        llr_v2[i] = llr2[i] + llr1[i] * v1_dec[i];
    }
    // Спускаемся вправо
    sc_recursive<T>(leaf_offset + half, llr_v2, v2_dec, u_bits, bit_mask, ws);

    // Схлопывание узла для передачи наверх родителю (x1 = v1 * v2, x2 = v2)
    for (size_t i = 0; i < half; ++i)
    {
        output_codeword[i] = v1_dec[i] * v2_dec[i];
        output_codeword[i + half] = v2_dec[i];
    }
}

/**
 * @brief Чистый верхнеуровневый интерфейс ванильного SC без постобработки
 */
template <std::floating_point T>
inline std::vector<int> polar_decode_sc_clean(std::span<const T> llr, const std::vector<uint8_t> &bit_mask, DecoderWorkspace<T> &ws)
{
    const size_t n = llr.size();
    std::vector<int> u(n); // Вектор u сразу под биты 0 и 1
    std::vector<T> dummy_codeword(n); // Временный буфер для корня дерева

    typename DecoderWorkspace<T>::Guard top_guard(ws);
    ws.current_offset = 0;

    // Запускаем — вся магия происходит внутри листьев рекурсии
    sc_recursive<T>(0, llr, dummy_codeword, u, bit_mask, ws);

    return u; // Возвращаем чистый, готовый вектор u!
}

/**
 * @brief Сжатие и извлечение полезной информационной последовательности по маске 5G.
 *
 * @param decoded_bits Полный вектор u длины N.
 * @param k Число полезных информационных битов K.
 * @param bit_mask Глобальная полярная маска.
 * @return std::vector<int> Чистый вектор информации длины K.
 */
inline std::vector<int> polar_extract_info_generic(const std::vector<int> &decoded_bits, int k, const std::vector<uint8_t> &bit_mask)
{
    std::vector<int> info;
    info.reserve(k);
    for (size_t i = 0; i < decoded_bits.size(); ++i)
    {
        if (bit_mask[i] == 1)
        {
            info.push_back(decoded_bits[i]);
        }
    }
    return info;
}
