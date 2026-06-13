import matplotlib
matplotlib.use('Agg')  # Отключает GUI для работы в WSL
import matplotlib.pyplot as plt

plt.rcParams['mathtext.fontset'] = 'cm'

fig, ax = plt.subplots(figsize=(7, 5), facecolor='white')

# Координаты узлов (Явные X и Y)
v_pos = {1: (1.5, 3.5), 2: (3.5, 3.5), 3: (5.5, 3.5)}
c_pos = {1: (1.5, 1.0), 2: (3.5, 1.0), 3: (5.5, 1.0)}

# 1. Отрисовка ребер цикла (выделяем красным цветом)
cycle_edges = [(1,1), (1,2), (2,1), (2,2)]
for v, c in cycle_edges:
    ax.plot([v_pos[v][0], c_pos[c][0]], [v_pos[v][1], c_pos[c][1]], 
            color='#d62728', linestyle='-', linewidth=2, zorder=1)

# Обычные ребра (нейтральный серый пунктир) - убрали дублирующий цикл (2,3)
other_edges = [(3,2), (3,3)]
for v, c in other_edges:
    ax.plot([v_pos[v][0], c_pos[c][0]], [v_pos[v][1], c_pos[c][1]], 
            color='gray', linestyle=':', linewidth=1.5, zorder=1)

# 2. Отрисовка битовых узлов (кружки)
for v, pos in v_pos.items():
    ax.scatter(pos[0], pos[1], s=700, color='#1f77b4', edgecolors='black', zorder=2)
    # Размер шрифта увеличен до 15 для лучшей читаемости
    ax.text(pos[0], pos[1], f'$v_{v}$', color='white', ha='center', va='center', 
            fontsize=15, fontweight='bold', zorder=3)

# 3. Отрисовка проверочных узлов (квадраты)
for c, pos in c_pos.items():
    ax.scatter(pos[0], pos[1], s=600, color='#9467bd', marker='s', edgecolors='black', zorder=2)
    # Размер шрифта увеличен до 15 для лучшей читаемости
    ax.text(pos[0], pos[1], f'$c_{c}$', color='white', ha='center', va='center', 
            fontsize=15, fontweight='bold', zorder=3)

# Настройки границ отображения
ax.set_xlim(0, 7)
ax.set_ylim(0, 5)
ax.axis('off')

# 4. Боковые заголовки слоев (верхний над узлами, нижний в самом низу под узлами)
ax.text(0.5, 4.2, 'Битовые узлы (LLR из канала):', fontsize=11, ha='left', va='center', 
        color='#1f77b4', fontweight='bold', zorder=3)
ax.text(0.5, 0.3, 'Проверочные узлы (Операция $\\boxplus$):', fontsize=11, ha='left', va='center', 
        color='#9467bd', fontweight='bold', zorder=3)

# Центральная выноска-пояснение для узелка зацикливания
ax.text(3.5, 2.25, 'Короткий цикл\nдлиной 4 (Girth = 4)\nЗацикливание LLR!', 
        color='#d62728', fontsize=10, ha='center', va='center', zorder=3,
        bbox=dict(boxstyle='round,pad=0.4', fc='#fbe7e7', ec='#d62728', alpha=0.95))

plt.savefig("tanner_graph_cycle.svg", bbox_inches='tight', facecolor='white')
print("Векторный файл 'tanner_graph_cycle.svg' успешно перезаписан с корректной геометрией.")
