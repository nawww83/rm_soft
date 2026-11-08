cmake_minimum_required(VERSION 3.20)
project(reed_muller_simulation CXX)

# Требуем стандарт C++20
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Включаем максимальную оптимизацию по скорости и предупреждения
if(MSVC)
    add_compile_options(/O2 /W4)
else()
    add_compile_options(-O3 -Wall -Wextra)
endif()

# Поиск библиотеки для многопоточности (std::async требует pthreads на Linux)
find_package(Threads REQUIRED)

# Добавляем исполняемый файл проекта
add_executable(rm_sim main.cpp rm_codes.hpp)

# Линкуем потоки
target_link_libraries(rm_sim PRIVATE Threads::Threads)
