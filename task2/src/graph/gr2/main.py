import matplotlib.pyplot as plt

# Данные
processes = [2, 4, 8, 16]
T1 = 56.9  # Время последовательной программы
Tp = [28.22, 25.9, 18.3, 9.9]  # Время параллельной программы для каждого количества процессов

# Рассчёт ускорения
speedup = [T1 / t for t in Tp]

# Построение графика
plt.plot(processes, speedup, 'o-', label='С использованием MPI', color='orange')

# Настройка осей и заголовков
plt.xlabel('Количество процессов')
plt.ylabel('Ускорение')
plt.title('Ускорение')
plt.legend()

# Включить сетку
plt.grid(True)

# Сохранение изображения графика
plt.savefig('speedup_graph.png')

# Можно также отобразить график, если нужно
# plt.show()

