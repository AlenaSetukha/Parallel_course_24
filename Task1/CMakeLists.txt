cmake_minimum_required(VERSION 3.16)# Проверка версии CMake.
									# Если версия установленной программы
									# старее указаной, произойдёт аварийный выход.

set(CMAKE_CXX_STANDARD 11)          # стандарт компилятора

project(SupComp_task1)			    # Название проекта


#compiler flags
set(CMAKE_CXX_FLAGS_DEBUG "-Og -g -Wall -Wextra -pedantic -fsanitize=leak, -fsanitize=undefinded, -fsanitize=memory")
set(CMAKE_CXX_FLAGS_RELEASE "-Ofast -flto -march=native")


# Find OpenMP package
find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

# Поиск пакета MPI
find_package(MPI REQUIRED)

include_directories(include)               # Где искать заголовочные файлы


set(SOURCE_LIB 
    src/coeff_calc.cpp
    src/geom_calc.cpp
    src/get_matrix.cpp
    src/Point.cpp
    src/write_toFile.cpp)
add_library(calc_lib STATIC ${SOURCE_LIB})	# Создание статической библиотеки с именем calc_lib


set(SOURCE_EXE main_MPI_2.cpp)		        # Установка переменной со списком исходников для исполняемого файла
add_executable(main_MPI ${SOURCE_EXE})		# Создает исполняемый файл с именем main_Sequent
target_link_libraries(main_MPI PUBLIC calc_lib)		# Линковка программы с библиотекой
# Связывание MPI с вашим исполняемым файлом
target_link_libraries(main_MPI PUBLIC MPI::MPI_CXX)




