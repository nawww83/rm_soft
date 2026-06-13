import matplotlib
matplotlib.use('Agg')  # Отключает GUI для работы в WSL

from matplotlib.ticker import LogLocator
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erfc

# Настройка шрифтов для корректного отображения математических символов в PDF/SVG
plt.rcParams['mathtext.fontset'] = 'cm'  # Computer Modern (классический стиль LaTeX)
plt.rcParams['svg.fonttype'] = 'none'   # Экспортировать текст как текст, а не как кривые

# ==========================================
# 1. Теоретические формулы
# ==========================================
def ber_theoretical_bpsk(ebn0_db):
    ebn0_linear = 10**(ebn0_db / 10.0)
    return 0.5 * erfc(np.sqrt(ebn0_linear))

def ber_theoretical_qam16(ebn0_db):
    ebn0_linear = 10**(ebn0_db / 10.0)
    def q_func(x):
        return 0.5 * erfc(x / np.sqrt(2.0))
    base_arg = np.sqrt(0.8 * ebn0_linear)
    return (3.0 / 4.0) * q_func(base_arg) + \
           (2.0 / 4.0) * q_func(3.0 * base_arg) - \
           (1.0 / 4.0) * q_func(5.0 * base_arg)

# ==========================================
# 2. Загрузка данных симулятора
# ==========================================
csv_filename_1 = "live_results_rm_2_6.csv"
try:
    data_rm = pd.read_csv(csv_filename_1)
except FileNotFoundError:
    print(f"Ошибка: Файл '{csv_filename_1}' не найден. Запустите симулятор.")
    exit()

csv_filename_2 = "live_results_polar_2_6.csv"
try:
    data_polar = pd.read_csv(csv_filename_2)
except FileNotFoundError:
    print(f"Ошибка: Файл '{csv_filename_2}' не найден. Запустите симулятор.")
    exit()

ebn0_sim_rm = data_rm["Eb/N0 (dB)"].values
ber_sim_rm = data_rm["BER"].values

mask_rm = ber_sim_rm > 0
ebn0_sim_rm_clean = ebn0_sim_rm[mask_rm]
ber_sim_rm_clean = ber_sim_rm[mask_rm]

ebn0_sim_polar = data_polar["Eb/N0 (dB)"].values
ber_sim_polar = data_polar["BER"].values

mask_polar = ber_sim_polar > 0
ebn0_sim_polar_clean = ebn0_sim_polar[mask_polar]
ber_sim_polar_clean = ber_sim_polar[mask_polar]

ebn0_dense = np.linspace(2.0, 12.0, 500)
ber_bpsk_theory = ber_theoretical_bpsk(ebn0_dense)
ber_qam16_theory = ber_theoretical_qam16(ebn0_dense)

# ==========================================
# 3. Отрисовка графиков (Оптимизировано для экранов)
# ==========================================
plt.figure(figsize=(9.5, 7), facecolor='white') # Немного увеличили высоту

# Теоретические кривые (пунктир)
plt.semilogy(ebn0_dense, ber_bpsk_theory, color='#2ca02c', linestyle='--', 
             label='Некодированная BPSK (Теория ФМн-2)', linewidth=1.5)
plt.semilogy(ebn0_dense, ber_qam16_theory, color='#9467bd', linestyle='--', 
             label='Некодированная 16-QAM (Теория 3GPP/Wi-Fi)', linewidth=1.5)

# Симуляционные данные
if len(ebn0_sim_rm_clean) > 0:
    plt.semilogy(ebn0_sim_rm_clean, ber_sim_rm_clean, color='#1f77b4', linestyle='-', marker='o',
                 label='Симуляция: RM(2, 6) + 16-QAM (Soft Recursive Decoding)', markersize=6, linewidth=2)

if len(ebn0_sim_polar_clean) > 0:
    plt.semilogy(ebn0_sim_polar_clean, ber_sim_polar_clean, color='#d62728', linestyle='-', marker='s',
                 label='Симуляция: Polar(64, 22) [3GPP 5G Mask] + Random Interleaver + 16-QAM (Fast-SSC)', markersize=6, linewidth=2)

plt.grid(True, which="both", linestyle=":", color='gray', alpha=0.5)

# Подписи осей
plt.xlabel('Отношение сигнал/шум $E_b$/$N_0$ (дБ)', fontsize=11)
plt.ylabel('Вероятность битовой ошибки BER', fontsize=11)
plt.title('Кривые помехоустойчивости цифровых систем связи', 
          fontsize=12, fontweight='bold', pad=15)

# === ОБНОВЛЕННЫЕ ГРАНИЦЫ ОСЕЙ ДЛЯ ОТОБРАЖЕНИЯ ТОЧКИ 11 дБ ===
plt.xlim(1.5, 12.5)
plt.ylim(1e-11, 1.0) # Опустили до 10^-11, чтобы точка 5.65*10^-10 была видна
plt.legend(fontsize=10, loc='lower left', framealpha=0.95)

# Сохраняем в векторном формате SVG
output_img = "ber_comparison.svg"


# Находим текущие оси графика
ax = plt.gca()

# Принудительно заставляем подписывать каждую степень 10 (base=10.0, subs=(1.0,), numticks=15)
ax.yaxis.set_major_locator(LogLocator(base=10.0, subs=(1.0,), numticks=15))

# Дополнительно можно включить мелкую (вспомогательную) сетку внутри каждого порядка
ax.yaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10)*0.1, numticks=15))
plt.grid(True, which="minor", linestyle=":", color='lightgray', alpha=0.3)


plt.savefig(output_img, bbox_inches='tight', facecolor='white')
print(f"Векторный график успешно обновлен: {output_img}")
