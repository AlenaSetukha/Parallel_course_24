import numpy as np
import matplotlib.pyplot as plt


def read_data(coord_file, value_file):
    # Считываем координаты
    with open(coord_file, 'r') as f:
        lines = f.readlines()
    N, M = map(int, lines[0].strip().split())
    coords = np.array([list(map(float, line.strip().split())) for line in lines[1:]])

    x = coords[:, 0]
    y = coords[:, 1]

    # Считываем значения функции
    with open(value_file, 'r') as f:
        lines = f.readlines()
    N_val, M_val = map(int, lines[0].strip().split())
    # Проверяем соответствие размеров
    if N != N_val or M != M_val:
        raise ValueError("Размеры N и M для координат и значений функции не совпадают.")
    
    values = np.array(list(map(float, lines[1:])))

    return x, y, values


def plot_function_on_grid(x, y, z):
    # Код функции остался без изменений (как в предыдущем ответе)
    if len(x) == 0 or len(y) == 0 or z.shape != (len(y), len(x)):
        raise ValueError("Неправильные размеры входных данных.")
    
    X, Y = np.meshgrid(x, y)

    plt.figure(figsize=(8, 6))
    contour = plt.contourf(X, Y, z, levels=20, cmap='viridis')
    plt.colorbar(contour, label='Значение функции')
    plt.title('Значения функции на двумерной сетке')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.grid(True)
    plt.show()

# Пример использования функций
coord_file = 'results/grid.txt'  # Файл с координатами (CSV)
value_file = 'results/res.txt'  # Файл со значениями функции (CSV)

x, y, z = read_data(coord_file, value_file)
plot_function_on_grid(x, y, z)