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
    unique_x = np.unique(x)
    unique_y = np.unique(y)
    
    X, Y = np.meshgrid(unique_x, unique_y)
    Z = z.reshape(len(unique_y), len(unique_x))

    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111, projection='3d')

    ax.plot_surface(X, Y, Z, cmap='viridis', edgecolor='none')

    # Настройка графика
    ax.set_title('Numerical Solution')
    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_zlabel('Z axis')

    # Показ графика
    plt.show()



# Пример использования функций
coord_file = 'results/grid.txt'  # Файл с координатами (CSV)
value_file = 'results/res.txt'  # Файл со значениями функции (CSV)
x, y, z = read_data(coord_file, value_file)
plot_function_on_grid(x, y, z)