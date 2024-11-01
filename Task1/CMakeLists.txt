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




include_directories(include)
set(SOURCE_LIB 
    src/coeff_calc.cpp
    src/geom_calc.cpp
    src/get_matrix.cpp
    src/Point.cpp
    src/write_toFile.cpp)			        # Тоже самое, но для библиотеки
add_library(calc_lib STATIC ${SOURCE_LIB})	# Создание статической библиотеки с именем calc_lib







set(SOURCE_EXE main_Sequent.cpp)		    # Установка переменной со списком исходников для исполняемого файла
add_executable(main_Sequent ${SOURCE_EXE})		# Создает исполняемый файл с именем main_Sequent
target_link_libraries(main_Sequent PUBLIC calc_lib)		# Линковка программы с библиотекой

#set(SOURCE_EXE main_OpenMP.cpp)		    # Установка переменной со списком исходников для исполняемого файла
#add_executable(main_OpenMP ${SOURCE_EXE})		# Создает исполняемый файл с именем main_OpenMP
#target_link_libraries(main_OpenMP PUBLIC calc_lib)		# Линковка программы с библиотекой