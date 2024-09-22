import matplotlib.pyplot as plt

# Данные из таблицы
chunk_sizes = ['дефолт', 100, 200]  # Размер чанка
static_times = [4.27951, 15.9614, 13.395]  # Время работы для static
dynamic_times = [135.495, 55.133, 52.8437]  # Время работы для dynamic
guided_times = [4.73732, 4.72655, 4.70682]  # Время работы для guided

# Построение графиков для каждого типа schedule
plt.plot(chunk_sizes, static_times, 'o-', label='static', color='gray')
plt.plot(chunk_sizes, dynamic_times, 'o-', label='dynamic', color='orange')
plt.plot(chunk_sizes, guided_times, 'o-', label='guided', color='blue')

# Настройка осей и заголовков
plt.xlabel('Размер чанка')
plt.ylabel('Время работы, сек')
plt.title('Зависимость времени от размера чанка')
plt.legend()

# Включить сетку
plt.grid(True)

# Сохранение графика
plt.savefig('chunk_time_graph.png')

# Можно также отобразить график
# plt.show()

