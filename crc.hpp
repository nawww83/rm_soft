#include <vector>
#include <cstdint>

class Crc8Engine {
public:
    // Полином CRC-8 ATM: x^8 + x^2 + x + 1 (0x07)
    static constexpr uint8_t POLYNOMIAL = 0x07;

    // Вычисляет 8-битную контрольную сумму для битового вектора
    static uint8_t calculate(const std::vector<int>& bits) {
        uint8_t crc = 0x00;
        for (int bit : bits) {
            uint8_t mix = (crc >> 7) ^ bit;
            crc <<= 1;
            if (mix & 1) {
                crc ^= POLYNOMIAL;
            }
        }
        return crc;
    }

    // Добавляет 8 бит CRC в конец информационного вектора
    static std::vector<int> append_crc(const std::vector<int>& info) {
        std::vector<int> tx_with_crc = info;
        uint8_t crc_val = calculate(info);
        
        // Разворачиваем байт CRC в 8 последовательных бит
        for (int i = 7; i >= 0; --i) {
            tx_with_crc.push_back((crc_val >> i) & 1);
        }
        return tx_with_crc;
    }

    // Проверяет, валиден ли вектор с CRC (остаток должен быть равен 0)
    static bool check(const std::vector<int>& bits_with_crc) {
        return calculate(bits_with_crc) == 0x00;
    }
};
