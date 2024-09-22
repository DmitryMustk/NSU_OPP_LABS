import matplotlib.pyplot as plt

# Данные с изображения
processes = [2, 4, 8, 16]
sequential_times = [56.9, 56.9, 56.9, 56.9]  # Время для последовательной программы фиксировано
parallel_times = [28.22, 25.9, 18.3, 9.9]  # Время для параллельной программы

# Построение графика
plt.plot(processes, sequential_times, 'o-', label='Последовательная программа', color='orange')
plt.plot(processes, parallel_times, 'o-', label='С использованием MPI', color='green')

# Настройка осей и заголовков
plt.xlabel('Количество процессов')
plt.ylabel('Время работы, сек')
plt.title('Зависимость по времени')
plt.legend()

# Включить сетку
plt.grid(True)

# Сохранение изображения в файл
plt.savefig('graph.png')

# Можно также отобразить график, если нужно
# plt.show()

