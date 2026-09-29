# rm_soft

LLR-based decoders for Reed-Muller and Polar codes in C++.

Реализация декодеров Reed-Muller и Polar кодов на C++ с использованием log-likelihood ratios (LLR).

## Описание

Этот проект представляет собой компактную C++-реализацию для моделирования и оценки декодеров кодов Рида–Маллера (Reed-Muller, RM) и Поляра (Polar) на основе LLR. Код рассчитан на исследовательские задачи, эксперименты и анализ характеристик BER/BLER в каналах с шумом.

## Что входит в проект

- генерация и декодирование кодов Рида–Маллера;
- реализация кодов Поляра и масок;
- мягкое декодирование на основе LLR;
- поддержка CRC;
- вспомогательные функции модуляции;
- инструменты для оценки ошибок при передаче;
- скрипты на Python для построения графиков BER/BLER.

## Основные особенности

- декодирование на основе soft decisions;
- поддержка RM и Polar code конструкций;
- анализ производительности в условиях зашумленного канала;
- удобство для исследования параметров кода и декодера;
- инструменты визуализации результатов.

## Структура проекта

```text
.
├── LICENSE
├── README.md
├── crc.hpp
├── generals.hpp
├── main.cpp
├── modulation.hpp
├── node_types.hpp
├── polar_decoder.hpp
├── polar_encoder.hpp
├── polar_mask.hpp
├── polar_scl_codec.hpp
├── rm_codes.hpp
├── nr_5g_polar_table.hpp
├── get_polar_mask_5g.py
├── plot_ber.py
├── plot_bler.py
├── tanner_cycle.py
└── ...
