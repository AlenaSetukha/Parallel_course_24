import numpy as np
import matplotlib.pyplot as plt

def plot_surface():
    # Генерация данных
    x = np.linspace(-1.5, 1.5, 100)
    y = np.linspace(-1.0, 1.0, 100)
    X, Y = np.meshgrid(x, y)

    # Вычисление Z в зависимости от условий
    Z = np.zeros_like(X)  # Изначально все значения Z равны нулю
    # Условие для значения Z внутри эллипса
    mask = X**2 + 4*Y**2 < 1
    Z[mask] = (1. - X[mask]**2 - 4. * Y[mask]**2) / 10.

    # Создание графика
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X, Y, Z, cmap='viridis', edgecolor='none')

    # Настройка графика
    ax.set_title('Analytical solution z = (1 - x^2 - 4y^2) / 10')
    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_zlabel('Z axis')

    # Показ графика
    plt.show()

# Вызов функции
plot_surface()