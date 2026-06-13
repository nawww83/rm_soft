#include <vector>
#include <span>
#include <algorithm>
#include <cmath>
#include <cassert>


template <typename T>
struct SclPath {
    T path_metric = 0.0;          // Метрика пути (чем меньше, тем вероятнее путь)
    std::vector<int> u_dec;       // Декодированные информационные и замороженные биты u
    bool is_active = false;       // Флаг активности траектории

    void reset(size_t n) {
        path_metric = 0.0;
        u_dec.assign(n, 0);
        is_active = false;
    }
};


template <typename T>
class PolarSclDecoder {
private:
    size_t N;
    size_t L; // Размер списка (например, 4 или 8)
    
    // Структуры для хранения промежуточных LLR на каждом уровне полярного дерева
    // Массив векторов LLR: уровень 'm' имеет размер N, уровень 'm-1' размер N/2 и т.д.
    std::vector<std::vector<std::vector<T>>> stage_llr;
    // Промежуточные жесткие решения (bit бит) для операции g
    std::vector<std::vector<std::vector<int>>> stage_bits;

    size_t m_stages;

public:
    PolarSclDecoder(size_t n, size_t list_size) : N(n), L(list_size) {
        m_stages = 0;
        while ((1ULL << m_stages) < N) m_stages++;

        // Аллокация памяти под L путей на каждом уровне декомпозиции Ариакана
        stage_llr.resize(L, std::vector<std::vector<T>>(m_stages + 1));
        stage_bits.resize(L, std::vector<std::vector<int>>(m_stages + 1));

        for (size_t l = 0; l < L; ++l) {
            for (size_t s = 0; s <= m_stages; ++s) {
                stage_llr[l][s].resize(1 << (m_stages - s), 0.0);
                stage_bits[l][s].resize(1 << (m_stages - s), 0);
            }
        }
    }

    // Главная функция декодирования
    std::vector<int> decode(std::span<const T> input_llr, const std::vector<uint8_t>& is_frozen) {
        std::vector<SclPath<T>> paths(L);
        paths[0].reset(N);
        paths[0].is_active = true;

        // Копируем входные LLR в нулевой уровень первой траектории
        std::copy(input_llr.begin(), input_llr.end(), stage_llr[0][0].begin());

        // Последовательно обходим все N подканалов (листьев дерева Ариакана)
        for (size_t step = 0; step < N; ++step) {
            // 1. Спускаемся по дереву вниз к текущему листу 'step', обновляя LLR (операция f)
            for (size_t l = 0; l < L; ++l) {
                if (paths[l].is_active) {
                    recursive_llr_update(l, step);
                }
            }

            // Проверяем, заморожен ли текущий подканал по маске 5G
            if (is_frozen[step] == 0) { 
                // --- ЗАМОРОЖЕННЫЙ КАНАЛ ---
                // Все активные пути обязаны принять решение '0'. Корректируем их метрики.
                for (size_t l = 0; l < L; ++l) {
                    if (paths[l].is_active) {
                        paths[l].u_dec[step] = 0;
                        T current_llr = stage_llr[l][m_stages][0];
                        // Если LLR отрицательный, а мы принудительно ставим 0, штрафуем метрику пути
                        if (current_llr < 0.0) {
                            paths[l].path_metric += std::abs(current_llr);
                        }
                        // Поднимаем жесткое решение по дереву вверх (операция g)
                        recursive_bits_update(l, step, 0);
                    }
                }
            } else {
                // --- ИНФОРМАЦИОННЫЙ КАНАЛ ---
                // Каждый активный путь раздваивается на две гипотезы: бит=0 и бит=1.
                // Всего получаем до 2*L возможных путей, из которых нужно выбрать ТОП-L лучших.
                std::vector<SclPath<T>> candidate_paths;
                candidate_paths.reserve(2 * L);

                for (size_t l = 0; l < L; ++l) {
                    if (paths[l].is_active) {
                        T current_llr = stage_llr[l][m_stages][0];

                        // Гипотеза 0: текущий бит равен 0
                        SclPath<T> p0 = paths[l];
                        p0.u_dec[step] = 0;
                        if (current_llr < 0.0) p0.path_metric += std::abs(current_llr);
                        candidate_paths.push_back(p0);

                        // Гипотеза 1: текущий бит равен 1
                        SclPath<T> p1 = paths[l];
                        p1.u_dec[step] = 1;
                        if (current_llr >= 0.0) p1.path_metric += std::abs(current_llr);
                        candidate_paths.push_back(p1);
                    }
                }

                // Сортируем кандидатов по возрастанию метрики пути (лучшие — в начале)
                std::sort(candidate_paths.begin(), candidate_paths.end(), [](const SclPath<T>& a, const SclPath<T>& b) {
                    return a.path_metric < b.path_metric;
                });

                // Оставляем только ТОП-L путей, остальные уничтожаем
                size_t active_count = std::min(L, candidate_paths.size());
                for (size_t i = 0; i < L; ++i) {
                    if (i < active_count) {
                        paths[i] = candidate_paths[i];
                        paths[i].is_active = true;
                        // Синхронизируем внутренние структуры LLR и бит для выживших путей
                        // (В полноценных кодеках тут копируют указатели, для симулятора проще скопировать массивы)
                        recursive_bits_update(i, step, paths[i].u_dec[step]);
                    } else {
                        paths[i].is_active = false;
                    }
                }
            }
        }

        // --- ВЕРИФИКАЦИЯ ЧЕРЕЗ CRC ---
        // Сканируем выжившие L путей от самых вероятных к наименее вероятным.
        // Первый же путь, который успешно проходит проверку контрольной суммы CRC, объявляется победителем!
        for (size_t l = 0; l < L; ++l) {
            if (paths[l].is_active) {
                if (Crc8Engine::check(paths[l].u_dec)) {
                    return paths[l].u_dec; // Нашли идеальный путь без ошибок!
                }
            }
        }

        // Если ни один путь не прошел CRC (тяжелая ошибка в канале), 
        // возвращаем просто самый первый, наиболее вероятный по метрике путь.
        return paths[0].u_dec;
    }

private:
    // Рекурсивный спуск по дереву Ариакана для пересчета LLR (Мягкие f и g операции)
    void recursive_llr_update(size_t l, size_t step) {
        size_t layer = 0;
        size_t temp_step = step;
        
        // Ищем уровень, на котором требуется обновление
        while (layer < m_stages && (temp_step % 2 == 0)) {
            size_t sub_block_size = 1 << (m_stages - layer);
            size_t half = sub_block_size / 2;
            size_t node_idx = temp_step / sub_block_size;
            
            auto src_llr = std::span<const T>(stage_llr[l][layer]).subspan(node_idx * sub_block_size, sub_block_size);
            auto dest_llr = std::span<T>(stage_llr[l][layer + 1]).subspan(node_idx * half, half);

            // Вычисляем операцию f для левой ветви
            for (size_t i = 0; i < half; ++i) {
                T l1 = src_llr[i];
                T l2 = src_llr[half + i];
                T s1 = (l1 >= 0.0) ? 1.0 : -1.0;
                T s2 = (l2 >= 0.0) ? 1.0 : -1.0;
                dest_llr[i] = s1 * s2 * std::min(std::abs(l1), std::abs(l2));
            }
            layer++;
            temp_step /= 2;
        }

        // Если мы находимся в правой ветви (индекс нечетный), выполняем операцию g
        if (layer < m_stages) {
            size_t sub_block_size = 1 << (m_stages - layer);
            size_t half = sub_block_size / 2;
            size_t node_idx = temp_step / sub_block_size;

            auto src_llr = std::span<const T>(stage_llr[l][layer]).subspan(node_idx * sub_block_size, sub_block_size);
            auto dest_llr = std::span<T>(stage_llr[l][layer + 1]).subspan(node_idx * half, half);
            auto bit_mask = std::span<const int>(stage_bits[l][layer + 1]).subspan(node_idx * half, half);

            for (size_t i = 0; i < half; ++i) {
                // Операция g: llr2 + llr1 * (1 - 2*bit)
                T bit_sign = (bit_mask[i] == 0) ? 1.0 : -1.0;
                dest_llr[i] = src_llr[half + i] + src_llr[i] * bit_sign;
            }
        }
    }

    // Подъем по дереву вверх для трансляции жестких решений
    void recursive_bits_update(size_t l, size_t step, int bit_val) {
        stage_bits[l][m_stages][0] = bit_val;
        size_t layer = m_stages;
        size_t temp_step = step;

        while (layer > 0 && (temp_step % 2 == 1)) {
            size_t sub_block_size = 1 << (m_stages - layer + 1);
            size_t half = sub_block_size / 2;
            size_t node_idx = temp_step / sub_block_size;

            auto child_bits = std::span<const int>(stage_bits[l][layer]).subspan(node_idx * sub_block_size, sub_block_size);
            auto parent_bits = std::span<int>(stage_bits[l][layer - 1]).subspan(node_idx * half, half);

            // Собираем бабочку Ариакана обратно: x1 = u1 ^ u2, x2 = u2
            for (size_t i = 0; i < half; ++i) {
                parent_bits[i] = child_bits[i] ^ child_bits[half + i];
                parent_bits[half + i] = child_bits[half + i];
            }
            layer--;
            temp_step /= 2;
        }
    }
};



// --- ВНУТРИ ТОЧКИ СИМУЛЯЦИИ ДЛЯ ПОЛЯРНОГО КОДА ---
int polar_k_total = get_rm_k(r, m); // Общий размер открытых каналов (например, 22)
int polar_k_info = polar_k_total - 8; // Полезные биты за вычетом CRC-8 (22 - 8 = 14)

std::vector<uint8_t> is_frozen = get_polar_5g_mask(n, polar_k_total);

// Инициализируем SCL-декодер со списком L=4 прямо перед циклом по пакетам
PolarSclDecoder<double> scl_decoder(n, 4); 

for (long long packet = 0; packet < max_packets; ++packet) {
    // Генерируем только полезные информационные биты (K - 8)
    std::vector<int> tx_info(polar_k_info);
    for (int& b : tx_info) b = bit_dist(gen);

    // Добавляем CRC-8 контрольную сумму (длина tx_info_with_crc становится равной K)
    std::vector<int> tx_info_with_crc = Crc8Engine::append_crc(tx_info);

    // Передаем в ваш рекурсивный полярный кодер (XOR на подъеме)
    std::vector<int> tx_codeword = pure_polar_encode(n, m, tx_info_with_crc, is_frozen);

    // ... модуляция, добавление шума, расчет LLR в padded_rx_llr ...

    // --- НА ПРИЕМЕ ---
    // Вызываем SCL декодер (он вернет вектор u длины N)
    std::vector<int> decoded_u = scl_decoder.decode(std::span<const double>(padded_rx_llr).subspan(0, n), is_frozen);

    // Экстрактор вытаскивает вектор u_info_with_crc (длиной K)
    std::vector<int> rx_info_with_crc = polar_extract_info_generic(decoded_u, polar_k_total, is_frozen);

    // Отрезаем CRC-8, оставляя только полезные информационные биты (первые K - 8 штук)
    std::vector<int> rx_info(rx_info_with_crc.begin(), rx_info_with_crc.begin() + polar_k_info);

    // Сравнение tx_info и rx_info на наличие ошибок, подсчет BER/BLER
}
