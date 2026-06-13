import matplotlib
matplotlib.use('Agg')  # Отключает GUI для работы в WSL

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erfc
from matplotlib.ticker import LogLocator

# Настройка шрифтов для корректного отображения в PDF/SVG
plt.rcParams['mathtext.fontset'] = 'cm'  # Computer Modern (классический стиль LaTeX)
plt.rcParams['svg.fonttype'] = 'none'   # Экспортировать текст как текст

# ==========================================
# 1. Теоретические формулы (BER -> BLER для N=64)
# ==========================================
N_block = 64

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

def bler_theoretical(ber_vector, n):
    return 1.0 - (1.0 - ber_vector)**n

# ==========================================
# 2. Загрузка данных симулятора
# ==========================================
csv_filename_1 = "live_results_rm_2_6.csv"
try:
    data_rm = pd.read_csv(csv_filename_1)
except FileNotFoundError:
    print(f"Ошибка: Файл '{csv_filename_1}' не найден.")
    exit()

csv_filename_2 = "live_results_polar_2_6.csv"
try:
    data_polar = pd.read_csv(csv_filename_2)
except FileNotFoundError:
    print(f"Ошибка: Файл '{csv_filename_2}' не найден.")
    exit()

# Данные Рида-Маллера
ebn0_sim_rm = data_rm["Eb/N0 (dB)"].values
bler_sim_rm = data_rm["BLER (FER)"].values

mask_rm = bler_sim_rm > 0
ebn0_sim_rm_clean = ebn0_sim_rm[mask_rm]
bler_sim_rm_clean = bler_sim_rm[mask_rm]

# Данные Полярного кода
ebn0_sim_polar = data_polar["Eb/N0 (dB)"].values
bler_sim_polar = data_polar["BLER (FER)"].values

mask_polar = bler_sim_polar > 0
ebn0_sim_polar_clean = ebn0_sim_polar[mask_polar]
bler_sim_polar_clean = bler_sim_polar[mask_polar]

# Расчет теории для блоков N=64
ebn0_dense = np.linspace(2.0, 12.0, 500)
bler_bpsk_theory = bler_theoretical(ber_theoretical_bpsk(ebn0_dense), N_block)
bler_qam16_theory = bler_theoretical(ber_theoretical_qam16(ebn0_dense), N_block)

# ==========================================
# 3. Отрисовка графиков BLER
# ==========================================
plt.figure(figsize=(9.5, 7), facecolor='white')

# Теоретические кривые BLER (пунктир)
plt.semilogy(ebn0_dense, bler_bpsk_theory, color='#2ca02c', linestyle='--', 
             label=f'Некодированная BPSK (Теория FER, N={N_block})', linewidth=1.5)
plt.semilogy(ebn0_dense, bler_qam16_theory, color='#9467bd', linestyle='--', 
             label=f'Некодированная 16-QAM (Теория FER, N={N_block})', linewidth=1.5)

# Симуляционные данные BLER
if len(ebn0_sim_rm_clean) > 0:
    plt.semilogy(ebn0_sim_rm_clean, bler_sim_rm_clean, color='#1f77b4', linestyle='-', marker='o',
                 label='Симуляция: RM(2, 6) + 16-QAM (Soft Recursive Decoding)', markersize=6, linewidth=2)

if len(ebn0_sim_polar_clean) > 0:
    plt.semilogy(ebn0_sim_polar_clean, bler_sim_polar_clean, color='#d62728', linestyle='-', marker='s',
                 label='Симуляция: Polar(64, 22) [5G Mask] + Random Interleaver + 16-QAM (Fast-SSC)', markersize=6, linewidth=2)

plt.grid(True, which="major", linestyle=":", color='gray', alpha=0.5)

# Оформление осей
plt.xlabel('Отношение сигнал/шум $E_b$/$N_0$ (дБ)', fontsize=11)
plt.ylabel('Вероятность блоковой ошибки BLER (FER)', fontsize=11)
plt.title('Кривые блоковой помехоустойчивости цифровых систем связи', 
          fontsize=12, fontweight='bold', pad=15)

plt.xlim(1.5, 12.5)
plt.ylim(1e-9, 1.5) # До 10^-9, так как минимальный ненулевой BLER полярного кода = 1.28*10^-8
plt.xticks(ebn0_sim_rm)

# --- ИСПРАВЛЕНИЕ СЕТКИ СТЕПЕНЕЙ 10 (ПОДПИСЬ КАЖДОГО ПОРЯДКА) ---
ax = plt.gca()
ax.yaxis.set_major_locator(LogLocator(base=10.0, subs=(1.0,), numticks=15))
ax.yaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10)*0.1, numticks=15))
plt.grid(True, which="minor", linestyle=":", color='lightgray', alpha=0.3)
# ---------------------------------------------------------------

plt.legend(fontsize=10, loc='lower left', framealpha=0.95)

# Сохранение в векторном формате SVG
output_img = "bler_comparison.svg"
plt.savefig(output_img, bbox_inches='tight', facecolor='white')
print(f"Векторный график BLER успешно создан: {output_img}")
