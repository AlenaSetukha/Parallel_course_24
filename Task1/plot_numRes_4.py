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


def plot_function_on_grid(x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4):

    unique_x1 = np.unique(x1)
    unique_y1 = np.unique(y1)

    X1, Y1 = np.meshgrid(unique_x1, unique_y1)
    Z1 = z1.reshape(len(unique_y1), len(unique_x1))


    unique_x2 = np.unique(x2)
    unique_y2 = np.unique(y2)

    X2, Y2 = np.meshgrid(unique_x2, unique_y2)
    Z2 = z2.reshape(len(unique_y2), len(unique_x2))

    unique_x3 = np.unique(x3)
    unique_y3 = np.unique(y3)

    X3, Y3 = np.meshgrid(unique_x3, unique_y3)
    Z3 = z3.reshape(len(unique_y3), len(unique_x3))

    unique_x4 = np.unique(x4)
    unique_y4 = np.unique(y4)

    X4, Y4 = np.meshgrid(unique_x4, unique_y4)
    Z4 = z4.reshape(len(unique_y4), len(unique_x4))


    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111, projection='3d')

    ax.plot_surface(X1, Y1, Z1, cmap='viridis', edgecolor='none')
    ax.plot_surface(X2, Y2, Z2, cmap='viridis', edgecolor='none')
    ax.plot_surface(X3, Y3, Z3, cmap='viridis', edgecolor='none')
    ax.plot_surface(X4, Y4, Z4, cmap='viridis', edgecolor='none')

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

    parser.add_argument('coord_file2', type=str, help='File containing coordinates (grid.txt)')
    parser.add_argument('value_file2', type=str, help='File containing function values (res.txt)')

    parser.add_argument('coord_file3', type=str, help='File containing coordinates (grid.txt)')
    parser.add_argument('value_file3', type=str, help='File containing function values (res.txt)')

    parser.add_argument('coord_file4', type=str, help='File containing coordinates (grid.txt)')
    parser.add_argument('value_file4', type=str, help='File containing function values (res.txt)')

    args = parser.parse_args()

    x1, y1, z1 = read_data(args.coord_file1, args.value_file1)
    x2, y2, z2 = read_data(args.coord_file2, args.value_file2)
    x3, y3, z3 = read_data(args.coord_file3, args.value_file3)
    x4, y4, z4 = read_data(args.coord_file4, args.value_file4)
    plot_function_on_grid(x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4)


if __name__ == "__main__":
    main()