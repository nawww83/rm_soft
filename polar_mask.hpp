#pragma once

#include <iostream>
#include <vector>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <span>
#include <iomanip>
#include <random>
#include <bitset>

struct PolarChannel {
    size_t index;       // Индекс
    double m_value;     // Математическое ожидание LLR
    double error_prob;  // Вероятность ошибки
};

class PolarDesignGA
{
private:
    // Аппроксимация Трифонова П.В.
    static inline double phi(double x)
    {
        if (x <= 0.0)   return 1.0;
        if (x > 35.0)   return 0.0; 

        if (x <= 10.0)
        {
            return std::exp(-0.4527 * std::pow(x, 0.86) + 0.0773);
        }
        else
        {
            return std::sqrt(M_PI / x) * std::exp(-x / 4.0) * (1.0 - 10.0 / (7.0 * x));
        }
    }

    static inline double phi_inv(double y)
    {
        if (y >= 0.999999) return 0.0;
        if (y <= 0.000001) return 35.0;

        double low = 0.0;
        double high = 35.0;
        
        for (int iter = 0; iter < 24; ++iter)
        {
            double mid = low + (high - low) / 2.0;
            if (phi(mid) < y)
            {
                high = mid; 
            }
            else
            {
                low = mid;
            }
        }
        return low;
    }

    // Рекурсивное ядро (Сверху-вниз / Natural Order)
    static void generate_recursive(std::span<PolarChannel> channels, double m_in) {
        const size_t n = channels.size();
        
        // Базовый случай: дошли до конкретного подканала
        if (n <= 1) {
            channels[0].m_value = m_in;
            return;
        }

        const size_t half = n / 2;

        // Расчет f-узла (ухудшение для левой ветви дерева)
        double m_f;
        if (m_in > 35.0) {
            m_f = m_in;
        } else {
            double phi_val = phi(m_in);
            m_f = phi_inv(1.0 - (1.0 - phi_val) * (1.0 - phi_val));
        }

        // Расчет g-узла (улучшение для правой ветви дерева)
        double m_g = 2.0 * m_in;
        if (m_g > 300.0) m_g = 300.0;

        // Спускаемся влево и вправо строго по топологии декодера
        generate_recursive(channels.subspan(0, half), m_f);
        generate_recursive(channels.subspan(half, half), m_g);
    }
public:
    static std::vector<PolarChannel> generate(size_t N, double snr_db) {
        double snr_linear = std::pow(10.0, snr_db / 10.0);
        double m_0 = 4.0 * snr_linear;

        std::vector<PolarChannel> channels(N);
        for (size_t i = 0; i < N; ++i) {
            channels[i].index = i;
            channels[i].m_value = 0.0;
            channels[i].error_prob = 0.0;
        }

        // Запуск рекурсивного построения "сверху-вниз"
        generate_recursive(channels, m_0);

        // Расчет финальных вероятностей ошибок по вычисленным m_value
        for (size_t i = 0; i < N; ++i) {
            if (channels[i].m_value <= 0.0) channels[i].error_prob = 0.5;
            else channels[i].error_prob = 0.5 * std::erfc(std::sqrt(channels[i].m_value) / 2.0);
        }

        return channels;
    }

    /**
    * Верификатор полярной последовательности.
    * Проверяет GA-расчет на частичный порядок.
    * @param channels - Вектор каналов.
    */
    static bool verify_sequence(const std::vector<PolarChannel> &channels)
    {
        const size_t N = channels.size();
        bool is_valid = true;
        const double EPSILON = 1e-5;

        // Частичный порядок Арикана (Arikan's Partial Ordering)
        for (size_t i = 0; i < N; ++i)
        {
            for (size_t j = i + 1; j < N; ++j)
            {
                if ((i & j) == i) // i является подмножеством j
                {
                    if (channels[i].m_value > channels[j].m_value + EPSILON)
                    {
                        std::cout << "[ОШИБКА ВЕРИФИКАЦИИ] Нарушен частичный порядок! "
                                << "Канал " << i << " (m=" << channels[i].m_value << ") "
                                << "надежнее канала " << j << " (m=" << channels[j].m_value << ")\n";
                        is_valid = false;
                    }
                }
            }
        }
        return is_valid;
    }

    /**
     * Вывод полярной последовательности в консоль
     * @param channels - Исходный вектор каналов от generate() (в естественном порядке)
     * @param is_info_mask - Маска информационных бит (размер N)
     * @param K - Количество информационных бит
     * @param snr_db - Расчетное SNR
     */
    static void log_sequence(const std::vector<PolarChannel>& channels, 
                             const std::vector<bool>& is_info_mask, 
                             size_t K, double snr_db) {
        const size_t N = channels.size();

        // 1. Шапка с параметрами кода
        std::cout << "========================================================================\n";
        std::cout << "                  ПОЛЯРНАЯ ПОСЛЕДОВАТЕЛЬНОСТЬ (МЕТОД GA)\n";
        std::cout << "========================================================================\n";
        std::cout << "Параметры: N = " << N << ", K = " << K 
                  << ", Скорость R = " << std::fixed << std::setprecision(4) << (static_cast<double>(K) / N)
                  << ", Design SNR = " << std::fixed << std::setprecision(2) << snr_db << " dB\n\n";

        // 2. Компактный вектор маски (0 - Frozen, 1 - Info) для быстрого копирования
        std::cout << "Маска подканалов в естественном порядке (0..N-1):\n";
        for (size_t i = 0; i < N; ++i) {
            std::cout << (is_info_mask[i] ? "1" : "0");
        }
        std::cout << "\n\n";

        // 3. Сортировка каналов по надежности (от худших к лучшим)
        // В теории полярных кодов принято выстраивать последовательность от самых слабых к сильным
        std::vector<PolarChannel> sorted_channels = channels;
        std::stable_sort(sorted_channels.begin(), sorted_channels.end(), [](const PolarChannel& a, const PolarChannel& b) {
            return a.m_value < b.m_value; 
        });

        // 4. Форматированный вывод таблицы
        // Определяем ширину битового представления в зависимости от N (для N=64 нужно 6 бит)
        size_t num_bits = std::log2(N);

        std::cout << std::left 
                  << std::setw(8)  << "Rank" 
                  << std::setw(10) << "Index" 
                  << std::setw(12) << "Binary" 
                  << std::setw(16) << "M[LLR]" 
                  << std::setw(20) << "Prob. error" 
                  << std::setw(10) << "Type" 
                  << "\n";
        std::cout << std::string(76, '-') << "\n";

        std::cout << std::fixed << std::setprecision(5);
        for (size_t rank = 0; rank < N; ++rank) {
            const auto& ch = sorted_channels[rank];
            std::string type_str = is_info_mask[ch.index] ? "INFO" : "FROZEN";
            
            // Перевод индекса в двоичную строку нужной длины
            std::string bin_str = std::bitset<64>(ch.index).to_string().substr(64 - num_bits);

            std::cout << std::left 
                      << std::setw(8)  << rank + 1
                      << std::setw(10) << ch.index 
                      << std::setw(12) << bin_str
                      << std::setw(16) << ch.m_value 
                      << std::setw(20) << ch.error_prob 
                      << std::setw(10) << type_str 
                      << "\n";
        }
        std::cout << std::string(76, '-') << "\n\n";
    }

    /**
    * Вывод текстового log-профиля надежности подканалов на основе вероятности ошибки
    * Метрика: -log10(error_prob) — чем длиннее полоса, тем надежнее канал
    * @param channels - Вектор каналов в естественном порядке (0..N-1)
    * @param is_info_mask - Маска информационных бит (размер N)
    */
    static void log_reliability_profile(const std::vector<PolarChannel>& channels, 
                                    const std::vector<bool>& is_info_mask) {
        const size_t N = channels.size();
        
        std::cout << "========================================================================\n";
        std::cout << "          LOG-ПРОФИЛЬ НАДЕЖНОСТИ ПОДКАНАЛОВ ПО ВЕРОЯТНОСТИ ОШИБКИ\n";
        std::cout << "========================================================================\n";
        std::cout << "Метрика шкалы: -log10(error_prob). Чем длиннее шкала, тем меньше ошибок.\n";
        std::cout << "Условные обозначения: [.] - FROZEN канал, [#] - INFO канал\n";
        std::string line_76(76, '-');
        std::cout << line_76 << "\n";

        // 1. Вычисляем значения -log10(Pe) для всех каналов и находим максимум
        std::vector<double> log_pe_values(N);
        double max_log_pe = 0.0;

        for (size_t i = 0; i < N; ++i) {
            double pe = channels[i].error_prob;
            
            if (pe < 1e-15) {
                pe = 1e-15; 
            }
            if (pe > 0.499) {
                pe = 0.499;
            }

            log_pe_values[i] = -std::log10(pe);
            if (log_pe_values[i] > max_log_pe) {
                max_log_pe = log_pe_values[i];
            }
        }

        // 2. Отрисовка гистограммы
        const size_t max_bar_width = 40; 
        size_t num_bits = static_cast<size_t>(std::log2(N));

        for (size_t i = 0; i < N; ++i) {
            // Извлекаем оригинальный физический индекс канала из структуры
            size_t actual_channel_idx = channels[i].index;
            // Перевод индекса в бинарную строку
            std::string bin_str = std::bitset<64>(actual_channel_idx).to_string().substr(64 - num_bits);
    
            size_t bar_length = 0;
            if (max_log_pe > 0.0) {
                bar_length = static_cast<size_t>((log_pe_values[i] / max_log_pe) * max_bar_width);
            }
            if (bar_length == 0) bar_length = 1;

            char symbol = is_info_mask[actual_channel_idx] ? '#' : '.';
            std::string bar(bar_length, symbol);

            std::cout << "Ch " << std::setw(3) << std::left << actual_channel_idx 
                    << " (" << bin_str << ") "
                    << "Pe: [" << std::scientific << std::setprecision(2) << channels[i].error_prob << "] "
                    << "| " << std::setw(max_bar_width) << std::left << bar 
                    << "\n";
        }
        std::cout << line_76 << "\n\n";
    }
};

class PolarEncoder {
public:
    static inline void polar_encode_recursive_core(std::span<int> codeword) {
        const size_t n = codeword.size();
        if (n <= 1) return;
        const size_t half = n / 2;
        for (size_t i = 0; i < half; ++i) {
            codeword[i] = codeword[i] ^ codeword[half + i];
        }
        polar_encode_recursive_core(codeword.subspan(0, half));
        polar_encode_recursive_core(codeword.subspan(half, half));
    }

    static std::vector<int> encode(const std::vector<int>& info_bits, const std::vector<bool>& is_info_mask) {
        const size_t N = is_info_mask.size();
        std::vector<int> codeword(N, 0);

        size_t info_idx = 0;
        for (size_t i = 0; i < N; ++i) {
            if (is_info_mask[i]) {
                codeword[i] = info_bits[info_idx++];
            }
        }
        polar_encode_recursive_core(codeword);
        return codeword;
    }
};

class PolarDecoderSC {
private:
    static inline double f_node(double lambda1, double lambda2) {
        double sign = (lambda1 < 0.0) ^ (lambda2 < 0.0) ? -1.0 : 1.0;
        return sign * std::min(std::abs(lambda1), std::abs(lambda2));
    }
    static inline double g_node(double lambda1, double lambda2, int u_pol) {
        return u_pol == 0 ? (lambda1 + lambda2) : (lambda2 - lambda1);
    }

    /**
     * @brief Рекурсивное ядро декодирования
     * @param output_codeword Сюда узел записывает свои закодированные промежуточные биты (x) для родителя
     */
    static void decode_recursive(std::span<const double> input_llr, 
                                 const std::vector<bool>& is_info,
                                 std::vector<int>& final_u,
                                 size_t global_offset,
                                 std::span<double> lr_buffer,
                                 std::span<int> x_buffer, // Буфер для промежуточных кодовых слов
                                 std::span<int> output_codeword) 
    {
        const size_t n = input_llr.size();  
        // Базовый случай: лист дерева решений
        if (n <= 1) {
            if (!is_info[global_offset]) {
                final_u[global_offset] = 0;
                output_codeword[0] = 0;
            } else {
                int decision = (input_llr[0] < 0.0) ? 1 : 0;
                final_u[global_offset] = decision;
                output_codeword[0] = decision;
            }
            return;
        }
        const size_t half = n / 2;

        // Нарезка предаллоцированной памяти
        auto llr_left   = lr_buffer.subspan(0, half);
        auto llr_right  = lr_buffer.subspan(half, half);
        auto x_left     = x_buffer.subspan(0, half);
        auto x_right    = x_buffer.subspan(half, half);

        // Сдвиг указателей арены памяти для более глубоких уровней рекурсии
        auto next_lr_buffer = lr_buffer.subspan(2 * half);
        auto next_x_buffer  = x_buffer.subspan(2 * half);

        // 1. Левый подблок: расчет f-узлов и спуск влево
        for (size_t i = 0; i < half; ++i) {
            llr_left[i] = f_node(input_llr[i], input_llr[half + i]);
        }
        decode_recursive(llr_left, is_info, final_u, global_offset, next_lr_buffer, next_x_buffer, x_left);

        // 2. Правый подблок: расчет g-узлов
        for (size_t i = 0; i < half; ++i) {
            llr_right[i] = g_node(input_llr[i], input_llr[half + i], x_left[i]);
        }
        decode_recursive(llr_right, is_info, final_u, global_offset + half, next_lr_buffer, next_x_buffer, x_right);

        // 3. Сборка вверх: родителю нужно отдать закодированное слово текущего узла
        // x1 = u1 ^ u2, x2 = u2
        for (size_t i = 0; i < half; ++i) {
            output_codeword[i] = x_left[i] ^ x_right[i];
            output_codeword[half + i] = x_right[i];
        }
    }

public:
    static std::vector<int> decode(const std::vector<double>& rx_llr, const std::vector<bool>& is_info) {
        const size_t N = rx_llr.size();
        std::vector<int> final_u(N, 0);
        // Выделение рабочей памяти (арены)
        std::vector<double> lr_workspace(2 * N, 0.0);
        std::vector<int> x_workspace(2 * N, 0);
        std::vector<int> root_output_codeword(N, 0);
        decode_recursive(rx_llr, is_info, final_u, 0, lr_workspace, x_workspace, root_output_codeword);
        return final_u;
    }
};

/*
int main() {
    const size_t N = 16;       
    const size_t K = 8;       
    const double snr_db_design = 2.0; 

    // 1. Генерируем каналы (внутри они рассчитываются послойно)
    std::vector<PolarChannel> channels = PolarDesignGA::generate(N, snr_db_design);
    PolarDesignGA::verify_sequence(channels);

    // 2. Создаем чистую маску размера N (по умолчанию все FROZEN = false)
    std::vector<bool> is_info_mask(N, false);

    // 3. Делаем глубокую копию для сортировки
    std::vector<PolarChannel> sorted_channels = channels;

    // 4. Сортируем строго по возрастанию m_value
    std::stable_sort(sorted_channels.begin(), sorted_channels.end(), [](const PolarChannel& a, const PolarChannel& b) {
        return a.m_value < b.m_value;
    });

    // 5. Выделяем K самых надежных каналов (они гарантированно в конце sorted_channels)
    for (size_t i = N - K; i < N; ++i) {
        // Извлекаем оригинальный физический индекс канала, который ему присвоил генератор
        size_t best_channel_index = sorted_channels[i].index; 
        is_info_mask[best_channel_index] = true;
    }

    // Выводим профиль и маску
    PolarDesignGA::log_sequence(channels, is_info_mask, K, snr_db_design);
    PolarDesignGA::log_reliability_profile(channels, is_info_mask);

    // Инициализация генераторов случайных чисел
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> bit_dist(0, 1);
    std::normal_distribution<double> gauss_dist(0.0, 1.0);

    // Параметры канала (фиксированные)
    double R = static_cast<double>(K) / N;
    const double snr_db = 4.0; 
    double snr_linear = std::pow(10.0, snr_db / 10.0);
    double sigma = std::sqrt(1.0 / (2.0 * R * snr_linear));

    std::cout << "========================================================================\n";
    std::cout << "          СТАРТ ЦИКЛА ПЕРЕДАЧИ (10 ИТЕРАЦИЙ, SNR = " << snr_db << " dB)\n";
    std::cout << "========================================================================\n\n";

    // 2. ЦИКЛ ПЕРЕДАЧИ КАДРОВ
    for (int iter = 1; iter <= 10; ++iter) {
        std::cout << "--- ИТЕРАЦИЯ №" << iter << " ---\n";

        // Генерация случайных информационных бит
        std::vector<int> tx_info(K);
        std::cout << "Исходные данные:     ";
        for (size_t i = 0; i < K; ++i) {
            tx_info[i] = bit_dist(gen);
            std::cout << tx_info[i];
        }
        std::cout << "\n";
        // Кодирование u -> x
        std::vector<int> codeword = PolarEncoder::encode(tx_info, is_info_mask);
        // Симуляция AWGN-канала (каждая итерация получит свою уникальную реализацию шума)
        std::vector<double> rx_llr(N);
        for (size_t i = 0; i < N; ++i) {
            double tx_symbol = (codeword[i] == 0) ? 1.0 : -1.0;
            double noise = gauss_dist(gen) * sigma;
            double rx_symbol = tx_symbol + noise;
            rx_llr[i] = 2.0 * rx_symbol / (sigma * sigma);
        }
        // Декодирование (SC возвращает напрямую вектор u)
        std::vector<int> decoded_u = PolarDecoderSC::decode(rx_llr, is_info_mask);
        // Извлечение инф. бит
        std::vector<int> rx_info;
        std::cout << "Декодированные данные: ";
        for (size_t i = 0; i < N; ++i) {
            if (is_info_mask[i]) {
                rx_info.push_back(decoded_u[i]);
                std::cout << decoded_u[i];
            }
        }
        std::cout << "\n";
        // Проверка на ошибки в текущем кадре
        size_t bit_errors = 0;
        for (size_t i = 0; i < K; ++i) {
            if (tx_info[i] != rx_info[i]) bit_errors++;
        }

        if (bit_errors == 0) {
            std::cout << "[УСПЕХ] Кадр №" << iter << " декодирован абсолютно верно без ошибок!\n";
        } else {
            std::cout << "[ОШИБКА] В кадре №" << iter << " обнаружено " << bit_errors << " битовых ошибок(и).\n";
        }
        std::cout << "-------------------------------------\n\n";
    }
    return 0;
}
*/