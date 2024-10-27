import numpy as np
import matplotlib.pyplot as plt

def plot_surface():
    # Генерация данных
    x = np.linspace(-1, 1, 100)
    y = np.linspace(-0.5, 0.5, 100)
    X, Y = np.meshgrid(x, y)
    Z = (1 - X**2 - 4*Y**2) / 10

    # Создание графика
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X, Y, Z, cmap='viridis', edgecolor='none')

    # Настройка графика
    ax.set_title('Numerical Solution')
    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_zlabel('Z axis')

    # Показ графика
    plt.show()

# Вызов функции
plot_surface()