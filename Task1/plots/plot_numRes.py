import numpy as np
import matplotlib.pyplot as plt
import argparse


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


def plot_function_on_grid(x1, y1, z1):

    unique_x1 = np.unique(x1)
    unique_y1 = np.unique(y1)

    X1, Y1 = np.meshgrid(unique_x1, unique_y1)
    Z1 = z1.reshape(len(unique_y1), len(unique_x1))


    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111, projection='3d')

    ax.plot_surface(X1, Y1, Z1, cmap='viridis', edgecolor='none')

    # Настройка графика
    ax.set_title('Numerical Solution')
    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_zlabel('Z axis')

    # Показ графика
    plt.show()


def main():
    parser = argparse.ArgumentParser(description='Plot function on grid from files.')   

    parser.add_argument('coord_file1', type=str, help='File containing coordinates (grid.txt)')
    parser.add_argument('value_file1', type=str, help='File containing function values (res.txt)')

    args = parser.parse_args()

    x1, y1, z1 = read_data(args.coord_file1, args.value_file1)
    plot_function_on_grid(x1, y1, z1)


if __name__ == "__main__":
    main()