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
        //  Функция phi(x) с защитой от переполнения.
        static inline double phi(double x)
        {
            if (x <= 0.0)
                return 1.0;
            if (x > 35.0)
                return 0.0; // Порог насыщения: при x > 35 значение phi(x) меньше предела точности double

            if (x <= 10.0)
            {
                return std::exp(-0.4527 * std::pow(x, 0.86) + 0.0218);
            }
            else
            {
                return std::sqrt(M_PI / x) * std::exp(-x / 4.0) * (1.0 - 10.0 / (7.0 * x));
            }
        }

        // Обратная функция phi^(-1)(y).
        static inline double phi_inv(double y)
        {
            if (y >= 0.9999999)
                return 0.0;
            if (y <= 0.0000001)
                return 35.0; // Если канал стал "идеальным", ограничиваем его рост порогом насыщения

            double low = 0.0;
            double high = 35.0; // Ищем в пределах разумного диапазона
            for (int iter = 0; iter < 64; ++iter)
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

        // Модифицированное ядро с защитой от вычислительного взрыва.
        static void ga_recursive_core_left_triangular(std::span<PolarChannel> channels)
        {
            const size_t n = channels.size();
            if (n <= 1)
                return;
            const size_t half = n / 2;

            for (size_t i = 0; i < half; ++i)
            {
                double m1 = channels[i].m_value;
                double m2 = channels[half + i].m_value;

                // f-узел: верхняя половина ухудшается за счет XOR
                double m_f;
                // Защита: если один из каналов уже "идеальный", операция XOR определяется вторым каналом
                if (m1 > 35.0)
                    m_f = m2;
                else if (m2 > 35.0)
                    m_f = m1;
                else
                {
                    m_f = phi_inv(1.0 - (1.0 - phi(m1)) * (1.0 - phi(m2)));
                }

                if (std::isnan(m_f) || std::isinf(m_f))
                    m_f = 0.0;

                // g-узел: нижняя половина улучшается
                double m_g = m1 + m2;

                // Ограничиваем рост LLR, чтобы избежать вычислительного хаоса на больших N
                if (m_g > 300.0)
                    m_g = 300.0;

                channels[i].m_value = m_f;
                channels[half + i].m_value = m_g;
            }

            ga_recursive_core_left_triangular(channels.subspan(0, half));
            ga_recursive_core_left_triangular(channels.subspan(half, half));
        }

    public:
        static std::vector<PolarChannel> generate(size_t N, double snr_db)
        {
            double snr_linear = std::pow(10.0, snr_db / 10.0);
            double m_0 = 4.0 * snr_linear;

            std::vector<PolarChannel> channels(N);
            for (size_t i = 0; i < N; ++i)
            {
                channels[i].index = i;
                channels[i].m_value = m_0;
                channels[i].error_prob = 0.0;
            }

            // Запуск расчета
            ga_recursive_core_left_triangular(channels);

            // Расчет вероятности ошибки (согласованный с m_0 = 4 * SNR)
            for (size_t i = 0; i < N; ++i)
            {
                if (channels[i].m_value <= 0.0)
                {
                    channels[i].error_prob = 0.5;
                }
                else
                {
                    channels[i].error_prob = 0.5 * std::erfc(std::sqrt(channels[i].m_value / 8.0));
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

            for (size_t i = 0; i < N; ++i)
            {
                for (size_t j = 0; j < N; ++j)
                {
                    // Проверяем отношение побитового включения: i является подмножеством j (i <= j по полупорядку)
                    // Пример: 1 (0001) является подмножеством 3 (0011). Значит канал 1 ОБЯЗАН быть хуже или равен каналу 3.
                    if ((i & j) == i)
                    {
                        if (channels[i].m_value > channels[j].m_value)
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
            std::sort(sorted_channels.begin(), sorted_channels.end(), [](const PolarChannel &a, const PolarChannel &b)
                      {
                          return a.m_value < b.m_value; // Сортировка по возрастанию надежности
                      });

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

            // 1. Сначала вычислим значения -log10(Pe) для всех каналов и найдем максимум
            std::vector<double> log_pe_values(N);
            double max_log_pe = 0.0;

            for (size_t i = 0; i < N; ++i)
            {
                double pe = channels[i].error_prob;

                // Защита от Pe = 0.0 (для идеальных каналов вроде 63-го при насыщении)
                if (pe < 1e-15)
                {
                    pe = 1e-15;
                }
                // Защита от худших каналов (если Pe близко к 0.5, log10 будет около -0.3)
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
            const size_t max_bar_width = 40; // Максимальная длина полосы
            size_t num_bits = std::log2(N);

            for (size_t i = 0; i < N; ++i)
            {
                // Перевод индекса в бинарную строку
                std::string bin_str = std::bitset<64>(channels[i].index).to_string().substr(64 - num_bits);

                // Пропорциональное вычисление длины полосы графика
                size_t bar_length = 0;
                if (max_log_pe > 0.0)
                {
                    bar_length = static_cast<size_t>((log_pe_values[i] / max_log_pe) * max_bar_width);
                }
                if (bar_length == 0)
                    bar_length = 1;

                char symbol = is_info_mask[i] ? '#' : '.';
                std::string bar(bar_length, symbol);

                std::cout << "Ch " << std::setw(3) << std::left << channels[i].index
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