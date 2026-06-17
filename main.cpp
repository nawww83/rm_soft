#include <iostream>
#include <vector>
#include <chrono>
#include <iomanip>
#include <cmath>
#include <random>
#include <complex>
#include <span>
#include <fstream>
#include <string>
#include <sstream>

#include "modulation.hpp"
#include "rm_codes.hpp"
#include "polar_encoder.hpp"
#include "polar_decoder.hpp"
#include "nr_5g_polar_table.hpp"

#include <future>
#include <thread>
#include <queue>
#include <mutex>
#include <condition_variable>
#include <functional>

// Простая и эффективная реализация Пула Потоков (Thread Pool)
class ThreadPool
{
private:
    std::vector<std::jthread> workers;
    std::queue<std::function<void()>> tasks;
    std::mutex queue_mutex;
    std::condition_variable cv;
    bool stop = false;

public:
    ThreadPool(size_t threads)
    {
        for (size_t i = 0; i < threads; ++i)
        {
            workers.emplace_back([this](std::stop_token st)
                                 {
                while (!st.stop_requested()) {
                    std::function<void()> task;
                    {
                        std::unique_lock<std::mutex> lock(this->queue_mutex);
                        this->cv.wait(lock, [this, &st] { 
                            return this->stop || !this->tasks.empty() || st.stop_requested(); 
                        });
                        
                        if ((this->stop || st.stop_requested()) && this->tasks.empty())
                            return;
                            
                        task = std::move(this->tasks.front());
                        this->tasks.pop();
                    }
                    task(); // Выполняем симуляцию точки
                } });
        }
    }

    // Шаблонная функция добавления задачи (возвращает std::future)
    template <class F, class... Args>
    auto enqueue(F &&f, Args &&...args)
        -> std::future<typename std::invoke_result<F, Args...>::type>
    {
        using return_type = typename std::invoke_result<F, Args...>::type;

        auto task = std::make_shared<std::packaged_task<return_type()>>(
            std::bind(std::forward<F>(f), std::forward<Args>(args)...));

        std::future<return_type> res = task->get_future();
        {
            std::lock_guard<std::mutex> lock(queue_mutex);
            if (stop)
                throw std::runtime_error("Добавление задачи в остановленный пул");
            tasks.emplace([task]()
                          { (*task)(); });
        }
        cv.notify_one();
        return res;
    }

    ~ThreadPool()
    {
        {
            std::lock_guard<std::mutex> lock(queue_mutex);
            stop = true;
        }
        cv.notify_all();
        // jthread автоматически сделают join при уничтожении вектора workers
    }
};

// Структура для возврата счетчиков ошибок за одну итерацию функции
struct SimResult
{
    double eb_n0_db;
    long long total_bits;
    long long error_bits;
    long long total_blocks;
    long long error_blocks;
};

// Структура для накопления всей статистики в main
struct TotalStats
{
    double eb_n0_db;
    long long total_bits = 0;
    long long error_bits = 0;
    long long total_blocks = 0;
    long long error_blocks = 0;
};

enum class CodeType
{
    ReedMuller,
    PolarFastSsc,
    PolarSc
};

class UniversalRandomInterleaver
{
private:
    size_t N;
    std::vector<size_t> p_forward;
    std::vector<size_t> p_backward;

public:
    explicit UniversalRandomInterleaver(size_t block_size, unsigned int seed = 42) : N(block_size)
    {
        p_forward.resize(N);
        std::iota(p_forward.begin(), p_forward.end(), 0);
        std::mt19937 g(seed);
        std::shuffle(p_forward.begin(), p_forward.end(), g);

        p_backward.resize(N);
        for (size_t i = 0; i < N; ++i)
        {
            p_backward[p_forward[i]] = i;
        }
    }

    size_t size() const { return N; }

    // ОПТИМИЗИРОВАНО: Принимает готовый вектор-приемник и перезаписывает его (без аллокаций!)
    inline void interleave(std::span<const int> src, std::span<int> dst) const
    {
        for (size_t i = 0; i < N; ++i)
        {
            dst[i] = src[p_forward[i]];
        }
    }

    // ОПТИМИЗИРОВАНО: Деперемежение LLR напрямую в существующий буфер
    template <std::floating_point T>
    inline void deinterleave(std::span<const T> src, std::span<T> dst) const
    {
        for (size_t i = 0; i < N; ++i)
        {
            dst[i] = src[p_backward[i]];
        }
    }
};

inline std::vector<uint8_t> get_polar_5g_mask(size_t n, int k) {
    assert(n <= 1024 && (n & (n - 1)) == 0 && "N должно быть степенью 2");
    std::vector<int> filtered;
    filtered.reserve(n);
    for (uint16_t idx : nr_5g_polar::reliability_table_1024) {
        if (idx < static_cast<int>(n)) {
            filtered.push_back(idx);
        }
    }
    std::vector<uint8_t> channel_mask(n, 0);
    size_t start_info_idx = filtered.size() - k;
    for (size_t i = start_info_idx; i < filtered.size(); ++i) {
        channel_mask[filtered[i]] = 1;
    }
    return channel_mask;
}

SimResult run_qam16_simulation_point(int r, int m, double eb_n0_db, CodeType code_type, int target_k = 0)
{
    const size_t n = 1 << m;

    // --- ОПРЕДЕЛЕНИЕ РАЗМЕРНОСТИ K И МАСКИ ЗАМОРОЗКИ ---
    int k = 0;
    std::vector<uint8_t> bit_mask;

    if (code_type == CodeType::ReedMuller)
    {
        k = get_rm_k(r, m);
    }
    else
    {
        k = (target_k == 0) ? get_rm_k(r, m) : target_k;
        // Получаем оригинальную маску 5G (1 = INFO, 0 = FROZEN)
        bit_mask = get_polar_5g_mask(n, k);
        // bit_mask = polar_mask::generate_custom_qam_mask(n, k, 5);
    }

    size_t padded_n = (n % 4 == 0) ? n : n + (4 - (n % 4));
    double code_rate = static_cast<double>(k) / n;
    double bits_per_symbol = 4.0; // 16-QAM
    double eb_n0_linear = std::pow(10.0, eb_n0_db / 10.0);

    double sigma2_total = 1.0 / (code_rate * bits_per_symbol * eb_n0_linear);
    double sigma2_1d = sigma2_total / 2.0;
    double sigma_1d = std::sqrt(sigma2_1d);

    auto hash_snr = std::hash<double>{}(eb_n0_db);
    auto hash_thread = std::hash<std::thread::id>{}(std::this_thread::get_id());
    auto time_now = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    std::mt19937 gen(static_cast<unsigned int>(hash_snr ^ hash_thread ^ time_now));

    std::uniform_int_distribution<int> bit_dist(0, 1);
    std::normal_distribution<double> noise_dist(0.0, sigma_1d);

    DecoderWorkspace<double> ws(n);
    UniversalRandomInterleaver interleaver(n, 42);

    long long total_bits = 0, error_bits = 0;
    long long total_blocks = 0, error_blocks = 0;

    const long long max_block_errors = (eb_n0_db >= 8.0) ? 2 : 20;
    const long long max_packets = (eb_n0_db >= 10.0) ? 300'000 : 50'000;

    std::vector<int> tx_info(k);
    std::vector<int> tx_codeword;
    std::vector<int> qam_tx_bits;
    qam_tx_bits.reserve(padded_n);
    std::vector<int> interleaved_bits_buf(n);     // Буфер для интерливера
    std::vector<double> deinterleaved_llr_buf(n); // Буфер для деинтерливера LLR

    for (long long packet = 0; packet < max_packets && error_blocks < max_block_errors; ++packet)
    {
        for (int &b : tx_info)
            b = bit_dist(gen);

        // --- ВЕТВЛЕНИЕ КОДИРОВАНИЯ И ПЕРЕМЕЖЕНИЯ ---
        if (code_type == CodeType::ReedMuller)
        {
            std::vector<int> tx_codeword = rm_encode(r, m, tx_info);
            qam_tx_bits.assign(padded_n, 0);
            std::copy(tx_codeword.begin(), tx_codeword.end(), qam_tx_bits.begin());
        }
        else if (code_type == CodeType::PolarFastSsc || code_type == CodeType::PolarSc)
        {
            std::vector<int> tx_codeword = pure_polar_encode(n, tx_info, bit_mask);
            interleaver.interleave(std::span<const int>(tx_codeword), std::span<int>(interleaved_bits_buf));
            qam_tx_bits.assign(padded_n, 0);
            std::copy(interleaved_bits_buf.begin(), interleaved_bits_buf.end(), qam_tx_bits.begin());
        }

        // --- МОДУЛЯЦИЯ (16-QAM) ---
        size_t num_symbols = padded_n / 4;
        std::vector<std::complex<double>> tx_symbols(num_symbols);
        for (size_t i = 0; i < num_symbols; ++i)
        {
            double i_ch = map_bit_pair_to_pam4(qam_tx_bits[i * 4 + 0], qam_tx_bits[i * 4 + 1]);
            double q_ch = map_bit_pair_to_pam4(qam_tx_bits[i * 4 + 2], qam_tx_bits[i * 4 + 3]);
            tx_symbols[i] = std::complex<double>(i_ch, q_ch);
        }

        // --- КАНАЛ С ШУМОМ ---
        std::vector<double> padded_rx_llr(padded_n);
        for (size_t i = 0; i < num_symbols; ++i)
        {
            double rx_i = tx_symbols[i].real() + noise_dist(gen);
            double rx_v = tx_symbols[i].imag() + noise_dist(gen);

            calculate_pam4_llr(rx_i, sigma2_1d, padded_rx_llr[i * 4 + 0], padded_rx_llr[i * 4 + 1]);
            calculate_pam4_llr(rx_v, sigma2_1d, padded_rx_llr[i * 4 + 2], padded_rx_llr[i * 4 + 3]);
        }

        // --- ВЕТВЛЕНИЕ ДЕПЕРЕМЕЖЕНИЯ И ДЕКОДИРОВАНИЯ ---
        std::vector<int> rx_info;
        if (code_type == CodeType::ReedMuller)
        {
            std::vector<int> rx_codeword = rm_soft_decode_fast(r, m, std::span<const double>(padded_rx_llr).subspan(0, n), ws);
            rx_info = rm_extract_info(r, m, rx_codeword);
        }
        else if (code_type == CodeType::PolarFastSsc)
        {
            std::span<const double> src_llr_span(padded_rx_llr.data(), n);
            std::span<double> dst_llr_span(deinterleaved_llr_buf.data(), n);
            interleaver.deinterleave<double>(src_llr_span, dst_llr_span);

            std::span<const double> ready_llr_span(deinterleaved_llr_buf.data(), n);
            std::vector<int> decoded_u = polar_decode_fast_generic<double>(ready_llr_span, bit_mask, ws);
            rx_info = polar_extract_info_generic(decoded_u, k, bit_mask);
        } else if (code_type == CodeType::PolarSc) {
             std::span<const double> src_llr_span(padded_rx_llr.data(), n);
            std::span<double> dst_llr_span(deinterleaved_llr_buf.data(), n);
            interleaver.deinterleave<double>(src_llr_span, dst_llr_span);

            std::span<const double> ready_llr_span(deinterleaved_llr_buf.data(), n);
            std::vector<int> decoded_u = polar_decode_sc_clean<double>(ready_llr_span, bit_mask, ws);
            rx_info = polar_extract_info_generic(decoded_u, k, bit_mask);
        }

        int errors_in_this_block = 0;
        for (int i = 0; i < k; ++i)
        {
            if (tx_info[i] != rx_info[i])
            {
                error_bits++;
                errors_in_this_block++;
            }
        }

        if (errors_in_this_block > 0)
        {
            error_blocks++;
        }

        total_bits += k;
        total_blocks++;
    }

    return {eb_n0_db, total_bits, error_bits, total_blocks, error_blocks};
}

int main()
{
    const CodeType code_type = CodeType::ReedMuller;
    const int r = 2;
    const int m = 6;

    {
        // =========================================================================
        // ЭКСПРЕСС-ТЕСТ (SELF-TEST) УНИВЕРСАЛЬНОГО ИНТЕРЛИВЕРА (IN-PLACE ВЕРСИЯ)
        // =========================================================================
        {
            const size_t test_n = 64; // Тестируем для N=64 (можно поменять на 128, 256 и т.д.)
            UniversalRandomInterleaver test_int(test_n, 42);

            // 1. Создаем тестовое кодовое слово из случайных бит 0 и 1
            std::vector<int> original_bits(test_n);
            std::mt19937 rand_gen(123);
            for (size_t i = 0; i < test_n; ++i)
            {
                original_bits[i] = rand_gen() % 2;
            }

            // 2. Предвыделяем буферы для работы (как в основном цикле симуляции)
            std::vector<int> interleaved_bits_buf(test_n);
            std::vector<double> deinterleaved_llr_buf(test_n);

            // 3. Выполняем перемешивание напрямую в буфер (in-place)
            test_int.interleave(original_bits, interleaved_bits_buf);

            // 4. Эмулируем прохождение канала (переводим перемешанные биты в LLR)
            std::vector<double> rx_llr_from_channel(test_n);
            for (size_t i = 0; i < test_n; ++i)
            {
                rx_llr_from_channel[i] = static_cast<double>(interleaved_bits_buf[i]);
            }

            // 5. Выполняем деперемежение напрямую в буфер (in-place)
            test_int.deinterleave<double>(rx_llr_from_channel, deinterleaved_llr_buf);

            // 6. Переводим восстановленные LLR обратно в жесткие биты для сверки
            std::vector<int> final_bits(test_n);
            for (size_t i = 0; i < test_n; ++i)
            {
                final_bits[i] = static_cast<int>(deinterleaved_llr_buf[i]);
            }

            // 7. Проверяем побитовое совпадение сквозного тракта
            bool test_success = true;
            for (size_t i = 0; i < test_n; ++i)
            {
                if (original_bits[i] != final_bits[i])
                {
                    std::cout << "ОШИБКА: Рассинхронизация in-place интерливера на индексе " << i
                              << " (Было: " << original_bits[i] << ", Стало: " << final_bits[i] << ")" << std::endl;
                    test_success = false;
                    break;
                }
            }

            if (test_success)
            {
                std::cout << "УСПЕХ: Оптимизированный двухтабличный интерливер для N=" << test_n
                          << " прошел верификацию бит-в-бит!" << std::endl;
            }
            else
            {
                std::cout << "Критическая ошибка интерливера!\n";
                return -1;
            }
        }
        // =========================================================================
    }

    std::stringstream ss;
    ss << "live_results_"
       << (code_type == CodeType::ReedMuller ? "rm" : (code_type == CodeType::PolarFastSsc ? "polar_fast_ssc" : "polar_sc"))
       << "_" << r << "_" << m << ".csv";

    const std::string csv_filename = ss.str();

    std::vector<double> ebn0_grid;
    for (double db = 2.0; db <= 12.0; db += 1.0)
    {
        ebn0_grid.push_back(db);
    }

    std::vector<TotalStats> accumulated_stats;
    for (double db : ebn0_grid)
    {
        accumulated_stats.push_back({db, 0, 0, 0, 0});
    }

    // Автоматически определяет количество логических ядер вашего CPU (например, 8, 12, 16)
    unsigned int num_threads = std::thread::hardware_concurrency();

    // Защита на случай, если функция вернула 0 (бывает на некоторых старых/специфичных ОС)
    if (num_threads == 0)
    {
        num_threads = 4;
    }

    std::cout << "[INFO] Обнаружено ядер CPU: " << num_threads << ". Инициализация пула потоков...\n";

    // Создаем пул потоков строго под архитектуру вашего процессора!
    ThreadPool pool(num_threads);

    std::cout << "Запуск непрерывной симуляции RM(" << r << ", " << m << ") + 16-QAM..." << std::endl;
    auto global_start_time = std::chrono::high_resolution_clock::now();
    long long current_iteration = 0;

    while (true)
    {
        current_iteration++;

        // 1. Асинхронный запуск расчетов по всем точкам SNR
        std::vector<std::future<SimResult>> futures;
        for (double db : ebn0_grid)
        {
            futures.push_back(pool.enqueue([=]()
                                           { return run_qam16_simulation_point(r, m, db, code_type); }));
        }
        // 2. Сбор результатов из потоков
        for (size_t i = 0; i < ebn0_grid.size(); ++i)
        {
            SimResult chunk = futures[i].get();
            accumulated_stats[i].total_bits += chunk.total_bits;
            accumulated_stats[i].error_bits += chunk.error_bits;
            accumulated_stats[i].total_blocks += chunk.total_blocks;
            accumulated_stats[i].error_blocks += chunk.error_blocks;
        }

        auto current_time = std::chrono::high_resolution_clock::now();
        auto elapsed_seconds = std::chrono::duration_cast<std::chrono::seconds>(current_time - global_start_time).count();

        // 3. Открываем файл для полной перезаписи актуальными данными
        std::ofstream csv_file(csv_filename);
        if (csv_file.is_open())
        {
            // Добавлена новая колонка "Simulation Time (s)" в заголовок CSV
            csv_file << "Eb/N0 (dB),BER,BLER (FER),Total Blocks,Error Blocks,Simulation Time (s)\n";
            csv_file << std::scientific << std::setprecision(6);
        }

        // 4. Очистка экрана консоли
        std::cout << "\033[2J\033[H";

        std::cout << "=================================================================\n";
        std::cout << "  УТОЧНЕНИЕ ХАРАКТЕРИСТИК ПОМЕХОУСТОЙЧИВОСТИ\n";
        std::cout << "=================================================================\n";
        std::cout << " Время работы: " << elapsed_seconds << " сек | Итерация обновления: " << current_iteration << "\n";
        std::cout << " Данные записываются в файл: " << csv_filename << "\n";
        std::cout << "-----------------------------------------------------------------\n";
        std::cout << " Eb/N0 (dB) |     BER    |    BLER (FER)  | Всего блоков (ошибок)\n";
        std::cout << "-----------------------------------------------------------------\n";

        // 5. Обработка, вывод на экран и запись в CSV
        for (const auto &stats : accumulated_stats)
        {
            double calculated_ber = (stats.total_bits > 0) ? static_cast<double>(stats.error_bits) / stats.total_bits : 0.0;

            double calculated_bler = (stats.total_blocks > 0) ? static_cast<double>(stats.error_blocks) / stats.total_blocks : 0.0;

            // Вывод в консоль
            std::cout << "    " << std::fixed << std::setprecision(1) << std::setw(4) << stats.eb_n0_db
                      << "    |  " << std::scientific << std::setprecision(4) << calculated_ber
                      << " |  " << std::scientific << std::setprecision(4) << calculated_bler
                      << " |  " << std::fixed << std::setprecision(0) << static_cast<double>(stats.total_blocks)
                      << " (" << std::fixed << std::setprecision(0) << stats.error_blocks << ")\n";

            // Запись текущей строки в CSV-файл (последним аргументом пишется elapsed_seconds)
            if (csv_file.is_open())
            {
                csv_file << std::fixed << std::setprecision(1) << stats.eb_n0_db << ","
                         << std::scientific << std::setprecision(6) << calculated_ber << ","
                         << calculated_bler << ","
                         << stats.total_blocks << ","
                         << stats.error_blocks << ","
                         << std::fixed << elapsed_seconds << "\n"; // Запись времени симуляции
            }
        }

        // Безопасное закрытие: только если файл был успешно открыт
        if (csv_file.is_open())
        {
            csv_file.close();
        }

        std::cout << "-----------------------------------------------------------------\n";
        std::cout << "Для прерывания симуляции нажмите Ctrl + C\n";

        std::this_thread::sleep_for(std::chrono::milliseconds(50));
    }

    return 0;
}
