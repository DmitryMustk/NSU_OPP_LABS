import matplotlib.pyplot as plt

# Данные
processes = [2, 4, 8, 16]
T1 = 56.9  # Время последовательной программы
Tp = [28.22, 25.9, 18.3, 9.9]  # Время параллельной программы для каждого количества процессов

# Рассчёт ускорения
speedup = [T1 / t for t in Tp]

# Рассчёт эффективности
efficiency = [(s / p) * 100 for s, p in zip(speedup, processes)]

# Построение графика
plt.plot(processes, efficiency, 'o-', label='Эффективность MPI', color='green')

# Настройка осей и заголовков
plt.xlabel('Количество процессов')
plt.ylabel('Эффективность (%)')
plt.title('Эффективность')
plt.legend()

# Включить сетку
plt.grid(True)

# Сохранение изображения графика
plt.savefig('efficiency_graph.png')

# Можно также отобразить график, если нужно
# plt.show()

