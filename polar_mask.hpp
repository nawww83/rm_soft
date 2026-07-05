#pragma once

#include <iostream>
#include <vector>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <span>
#include <iomanip>
#include <bitset>
#include <sstream>

namespace polar_mask
{
    // Структура подканала
    struct PolarChannel
    {
        size_t index;      // Индекс
        double m_value;    // Математическое ожидание LLR
        double error_prob; // Вероятность ошибки
    };

    class PolarDesignGA
    {
    private:
        /**
         * Аппроксимация Трифонова П.В. Адаптация метода Чунга (Chung) к полярным кодам.
         */
        static inline double phi(double x)
        {
            if (x <= 0.0)
                return 1.0;
            if (x > 35.0)
                return 0.0;

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
            if (y >= 0.999999)
                return 0.0;
            if (y <= 0.000001)
                return 35.0;

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

        // Расчет по слоям графа
        static void calculate_ga_stages(std::vector<PolarChannel> &channels, double m_0)
        {
            const size_t N = channels.size();

            // На нулевом этапе (входе) все каналы инициализируются начальным m_0
            for (size_t i = 0; i < N; ++i)
            {
                channels[i].m_value = m_0;
            }

            // Внешний цикл: шагаем по этапам (каскадам) графа Арикана.
            // Всего log2(N) этапов. step принимает значения: 1, 2, 4, 8 ... до N/2
            for (size_t step = 1; step < N; step <<= 1)
            {
                // Разделяем массив на независимые блоки размером 2 * step
                for (size_t block = 0; block < N; block += 2 * step)
                {
                    // Итерируемся внутри одного блока
                    for (size_t i = 0; i < step; ++i)
                    {
                        // Индексы гарантированно лежат в пределах [0, N-1]
                        size_t idx_left = block + i;
                        size_t idx_right = block + i + step;

                        double m_prev = channels[idx_left].m_value; // В начале шага они равны

                        // f-узел (левая ветвь) -> ухудшение
                        double m_f;
                        if (m_prev > 35.0)
                        {
                            m_f = m_prev;
                        }
                        else
                        {
                            double phi_val = phi(m_prev);
                            m_f = phi_inv(1.0 - (1.0 - phi_val) * (1.0 - phi_val));
                        }

                        // g-узел (правая ветвь) -> удвоение
                        double m_g = 2.0 * m_prev;
                        if (m_g > 300.0)
                            m_g = 300.0;

                        // Записываем результаты на свои места в текущем слое
                        channels[idx_left].m_value = m_f;
                        channels[idx_right].m_value = m_g;
                    }
                }
            }
        }

    public:
        static std::vector<PolarChannel> generate(size_t N, double snr_db)
        {
            // Проверка на степень двойки для защиты от бесконечных циклов
            if (N == 0 || (N & (N - 1)) != 0)
            {
                throw std::invalid_argument("Размер N должен быть степенью двойки.");
            }

            double snr_linear = std::pow(10.0, snr_db / 10.0);
            double m_0 = 4.0 * snr_linear;

            std::vector<PolarChannel> channels(N);
            for (size_t i = 0; i < N; ++i)
            {
                channels[i].index = i;
                channels[i].m_value = 0.0;
                channels[i].error_prob = 0.0;
            }

            // Запуск расчета
            calculate_ga_stages(channels, m_0);

            // Расчет финальной вероятности ошибки
            for (size_t i = 0; i < N; ++i)
            {
                if (channels[i].m_value <= 0.0)
                {
                    channels[i].error_prob = 0.5;
                }
                else
                {
                    channels[i].error_prob = 0.5 * std::erfc(std::sqrt(channels[i].m_value) / 2.0);
                }
            }

            return channels;
        }

        /**
         * Верификатор полярной последовательности.
         * Проверяет GA-расчет на соответствие фундаментальному свойству частичного порядка (Partial Ordering).
         * @param channels - Вектор каналов в естественном порядке индексов (0..N-1)
         * @return bool - true, если последовательность математически корректна
         */
        static bool verify_sequence(const std::vector<PolarChannel> &channels)
        {
            const size_t N = channels.size();
            bool is_valid = true;

            // Допуск для сравнения double (защита от погрешностей округления в зоне насыщения)
            const double EPSILON = 1e-6;

            // Оптимизация: j всегда строго больше i, так как подмножество по индексам
            // может приводить к большему числу только при добавлении единичных битов (i < j)
            for (size_t i = 0; i < N; ++i)
            {
                for (size_t j = i + 1; j < N; ++j)
                {
                    // Проверяем отношение побитового включения: i является подмножеством j
                    if ((i & j) == i)
                    {
                        // Использован допуск EPSILON, чтобы избежать ложных срабатываний
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
        static void log_sequence(const std::vector<PolarChannel> &channels,
                                 const std::vector<bool> &is_info_mask,
                                 size_t K, double snr_db)
        {
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
            for (size_t i = 0; i < N; ++i)
            {
                std::cout << (is_info_mask[i] ? "1" : "0");
            }
            std::cout << "\n\n";

            // 3. Сортировка каналов по надежности (от худших к лучшим)
            // В теории полярных кодов принято выстраивать последовательность от самых слабых к сильным
            std::vector<PolarChannel> sorted_channels = channels;
            std::stable_sort(sorted_channels.begin(), sorted_channels.end(), [](const PolarChannel &a, const PolarChannel &b)
                             { return a.m_value < b.m_value; });

            // 4. Форматированный вывод таблицы
            // Определяем ширину битового представления в зависимости от N (для N=64 нужно 6 бит)
            size_t num_bits = std::log2(N);

            std::cout << std::left
                      << std::setw(8) << "Rank"
                      << std::setw(10) << "Index"
                      << std::setw(12) << "Binary"
                      << std::setw(16) << "M[LLR]"
                      << std::setw(20) << "Prob. error"
                      << std::setw(10) << "Type"
                      << "\n";
            std::cout << std::string(76, '-') << "\n";

            std::cout << std::fixed << std::setprecision(5);
            for (size_t rank = 0; rank < N; ++rank)
            {
                const auto &ch = sorted_channels[rank];
                std::string type_str = is_info_mask[ch.index] ? "INFO" : "FROZEN";

                // Перевод индекса в двоичную строку нужной длины
                std::string bin_str = std::bitset<64>(ch.index).to_string().substr(64 - num_bits);

                std::cout << std::left
                          << std::setw(8) << rank + 1
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
        static void log_reliability_profile(const std::vector<PolarChannel> &channels,
                                            const std::vector<bool> &is_info_mask)
        {
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

            for (size_t i = 0; i < N; ++i)
            {
                double pe = channels[i].error_prob;

                if (pe < 1e-15)
                {
                    pe = 1e-15;
                }
                if (pe > 0.499)
                {
                    pe = 0.499;
                }

                log_pe_values[i] = -std::log10(pe);
                if (log_pe_values[i] > max_log_pe)
                {
                    max_log_pe = log_pe_values[i];
                }
            }

            // 2. Отрисовка гистограммы
            const size_t max_bar_width = 40;
            size_t num_bits = static_cast<size_t>(std::log2(N));

            for (size_t i = 0; i < N; ++i)
            {
                // Извлекаем оригинальный физический индекс канала из структуры
                size_t actual_channel_idx = channels[i].index;

                // Перевод оригинального индекса в бинарную строку
                std::string bin_str = std::bitset<64>(actual_channel_idx).to_string().substr(64 - num_bits);

                size_t bar_length = 0;
                if (max_log_pe > 0.0)
                {
                    bar_length = static_cast<size_t>((log_pe_values[i] / max_log_pe) * max_bar_width);
                }
                if (bar_length == 0)
                    bar_length = 1;

                // Берем статус "заморозки" строго по оригинальному индексу канала
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

} // namespace polar_mask

/*
int main()
{
    using namespace polar_mask;

    const size_t N = 16;
    const size_t K = 8;
    const double snr_db = 2.0;

    // Генерируем каналы (получаем вектор из N элементов)
    std::vector<PolarChannel> channels = PolarDesignGA::generate(N, snr_db);

    // Делаем копию вектора для сортировки, чтобы не нарушить порядок оригинального вектора
    std::vector<PolarChannel> sorted = channels;

    // Сортируем каналы по убыванию надежности (самые надежные с максимальным m_value — в начале)
    std::sort(sorted.begin(), sorted.end(), [](const PolarChannel &a, const PolarChannel &b)
              { return a.m_value > b.m_value; });

    // Создаем маску из N элементов, изначально заполненную false (все FROZEN)
    std::vector<bool> is_info_mask(N, false);

    // Берем первые K самых надежных каналов из отсортированного списка и помечаем их как true (INFO)
    for (size_t i = 0; i < K; ++i)
    {
        size_t best_index = sorted[i].index; // Узнаем исходный индекс этого канала
        is_info_mask[best_index] = true;     // Выделяем его под информационный бит
    }

    PolarDesignGA::log_sequence(channels, is_info_mask, K, snr_db);
    PolarDesignGA::log_reliability_profile(channels, is_info_mask);

    return 0;
}
*/