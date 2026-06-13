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
inline void pure_polar_decode_recursive(size_t leaf_offset, std::span<const T> llr, 
                                       std::span<T> output_codeword, 
                                       const std::vector<uint8_t>& bit_mask, DecoderWorkspace<T>& ws) 
{
    const size_t n = llr.size();
    NodeType type = determine_node_type(leaf_offset, n, bit_mask);
    
    // =========================================================================
    // БЫСТРЫЕ УЗЛЫ ОСТАНОВА (Fast-SSC Оптимизации)
    // =========================================================================
    if (type == NodeType::Rate0) {
        // Все биты поддерева заморожены (u = 0). Согласно 3GPP, логический 0 отображается в +1.0
        std::fill(output_codeword.begin(), output_codeword.end(), T(1)); 
        return; 
    }
    if (type == NodeType::Rate1) {
        // Узел чистой информации. Выполняется жесткое решение по знаку LLR
        for (size_t i = 0; i < n; ++i) {
            output_codeword[i] = (llr[i] >= T(0)) ? T(1) : T(-1); 
        }
        return;
    }
    if (type == NodeType::Rep) {
        // Узел повторения. Решение принимается по знаку суммы всех LLR блока
        T sum_llr = 0;
        for (size_t i = 0; i < n; ++i) sum_llr += llr[i];
        T bit_sign = (sum_llr >= T(0)) ? T(1) : T(-1);
        std::fill(output_codeword.begin(), output_codeword.end(), bit_sign); 
        return;
    }
    if (n == 1) {
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
    for (size_t i = 0; i < half; ++i) {
        T s1 = (llr1[i] >= T(0)) ? T(1) : T(-1);
        T s2 = (llr2[i] >= T(0)) ? T(1) : T(-1);
        llr_v1[i] = s1 * s2 * std::min(std::abs(llr1[i]), std::abs(llr2[i]));
    }
    pure_polar_decode_recursive<T>(leaf_offset, llr_v1, v1_dec, bit_mask, ws);
    
    auto llr_v2 = ws.allocate(half);
    auto v2_dec = ws.allocate(half); 
    
    // Правая ветвь: g-операция 3GPP (знак ПЛЮС согласован со спецификацией TS 38.212)
    for (size_t i = 0; i < half; ++i) {
        llr_v2[i] = llr2[i] + llr1[i] * v1_dec[i];
    }
    pure_polar_decode_recursive<T>(leaf_offset + half, llr_v2, v2_dec, bit_mask, ws);
    
    // Объединение решений (схлопывание узла Плоткина: x1 = v1 * v2, x2 = v2)
    for (size_t i = 0; i < half; ++i) {
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
inline std::vector<int> polar_decode_fast_generic(std::span<const T> llr, const std::vector<uint8_t>& bit_mask, DecoderWorkspace<T>& ws) {
    const size_t n = llr.size();
    std::vector<T> bipolar_output(n);
    
    typename DecoderWorkspace<T>::Guard top_guard(ws);
    ws.current_offset = 0; 
    
    pure_polar_decode_recursive<T>(0, llr, bipolar_output, bit_mask, ws);
    
    // Перевод из биполярного вида (+1/-1) в логические биты 5G (0/1)
    std::vector<int> codeword(n);
    for (size_t i = 0; i < n; ++i) {
        codeword[i] = (bipolar_output[i] > 0) ? 0 : 1;
    }
    
    // Прямой пошаговый разбор кодового слова x на исходные компоненты u (разбиение Плоткина).
    // Гарантирует идеальное сохранение индексов 5G маски.
    std::vector<int> u(n);
    auto extract_u_direct = [&](auto& self, size_t offset, std::span<const int> cw, std::span<int> out_u) -> void {
        size_t len = cw.size();
        if (len == 1) {
            out_u[offset] = cw[0];
            return;
        }
        size_t h = len / 2;
        std::vector<int> u1(h), u2(h);
        for(size_t i = 0; i < h; ++i) {
            u2[i] = cw[h + i];
            u1[i] = cw[i] ^ cw[h + i];
        }
        self(self, offset, u1, out_u);
        self(self, offset + h, u2, out_u);
    };
    
    extract_u_direct(extract_u_direct, 0, codeword, u);
    return u;
}

/**
 * @brief Сжатие и извлечение полезной информационной последовательности по маске 5G.
 * 
 * @param decoded_bits Полный вектор u длины N.
 * @param k Число полезных информационных битов K.
 * @param bit_mask Глобальная полярная маска.
 * @return std::vector<int> Чистый вектор информации длины K.
 */
inline std::vector<int> polar_extract_info_generic(const std::vector<int>& decoded_bits, int k, const std::vector<uint8_t>& bit_mask) {
    std::vector<int> info;
    info.reserve(k); 
    for (size_t i = 0; i < decoded_bits.size(); ++i) {
        if (bit_mask[i] == 1) {
            info.push_back(decoded_bits[i]);
        }
    }
    return info;
}
