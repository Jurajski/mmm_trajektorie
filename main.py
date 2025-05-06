import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import tkinter as tk

# Funkcje generujące różne sygnały wejściowe
def rectangular_signal(t, duration, amplitude):
    return amplitude if 0 <= t <= duration else 0

def triangular_signal(t, duration, amplitude):
    if 0 <= t <= duration:
        return (2 * amplitude / duration) * (duration/2 - abs(t - duration/2))
    return 0

def harmonic_signal(t, frequency, amplitude):
    return amplitude * np.sin(2 * np.pi * frequency * t)

# Funkcja wymuszająca
def F(t, signal_type, params):
    if signal_type == 'prostokątny':
        return rectangular_signal(t, params['duration'], params['amplitude'])
    elif signal_type == 'trójkątny':
        return triangular_signal(t, params['duration'], params['amplitude'])
    elif signal_type == 'harmoniczny':
        return harmonic_signal(t, params['frequency'], params['amplitude'])

# Metoda Eulera
def euler_method(x0, v0, t0, tf, dt, signal_type, params, M, b, k):
    n = int((tf - t0) / dt)
    t = np.linspace(t0, tf, n)
    x = np.zeros(n)
    v = np.zeros(n)
    x[0] = x0
    v[0] = v0

    for i in range(1, n):
        a = (F(t[i-1], signal_type, params) - b * v[i-1] - k * x[i-1]) / M
        v[i] = v[i-1] + a * dt
        x[i] = x[i-1] + v[i-1] * dt

    return t, x, v

# Metoda Rungego-Kutty 4-go rzędu (RK4)
def rk4_method(x0, v0, t0, tf, dt, signal_type, params, M, b, k):
    n = int((tf - t0) / dt)
    t = np.linspace(t0, tf, n)
    x = np.zeros(n)
    v = np.zeros(n)
    x[0] = x0
    v[0] = v0

    for i in range(1, n):
        k1_v = (F(t[i-1], signal_type, params) - b * v[i-1] - k * x[i-1]) / M
        k1_x = v[i-1]

        k2_v = (F(t[i-1] + dt/2, signal_type, params) - b * (v[i-1] + k1_v * dt/2) - k * (x[i-1] + k1_x * dt/2)) / M
        k2_x = v[i-1] + k1_v * dt/2

        k3_v = (F(t[i-1] + dt/2, signal_type, params) - b * (v[i-1] + k2_v * dt/2) - k * (x[i-1] + k2_x * dt/2)) / M
        k3_x = v[i-1] + k2_v * dt/2

        k4_v = (F(t[i-1] + dt, signal_type, params) - b * (v[i-1] + k3_v * dt) - k * (x[i-1] + k3_x * dt)) / M
        k4_x = v[i-1] + k3_v * dt

        v[i] = v[i-1] + (k1_v + 2 * k2_v + 2 * k3_v + k4_v) * dt / 6
        x[i] = x[i-1] + (k1_x + 2 * k2_x + 2 * k3_x + k4_x) * dt / 6

    return t, x, v

# Funkcja do aktualizacji wykresu
def update_plot(F_const, M, b, k, signal_type, duration, frequency):
    params = {'amplitude': F_const, 'duration': duration, 'frequency': frequency}
    t0, tf, dt = 0.0, 5.0, 0.05  # Krótsza symulacja z większym dt
    x0, v0 = 0.0, 0.0

    # Symulacja metodą Eulera
    t_euler, x_euler, v_euler = euler_method(x0, v0, t0, tf, dt, signal_type, params, M, b, k)
    
    # Symulacja metodą RK4
    t_rk4, x_rk4, v_rk4 = rk4_method(x0, v0, t0, tf, dt, signal_type, params, M, b, k)

    # Aktualizacja wykresów
    ax1.clear()
    ax1.plot(t_euler, v_euler, label="Euler", linestyle='--')
    ax1.plot(t_rk4, v_rk4, label="RK4", linestyle='-')
    ax1.set_title("Prędkość wózka")
    ax1.set_xlabel("Czas [s]")
    ax1.set_ylabel("Prędkość [m/s]")
    ax1.legend()
    ax1.grid(True)

    ax2.clear()
    ax2.plot(t_euler, x_euler, label="Euler", linestyle='--')
    ax2.plot(t_rk4, x_rk4, label="RK4", linestyle='-')
    ax2.set_title("Położenie wózka")
    ax2.set_xlabel("Czas [s]")
    ax2.set_ylabel("Położenie [m]")
    ax2.legend()
    ax2.grid(True)

    canvas.draw()

# Funkcja zamykająca aplikację
def on_closing():
    plt.close('all')
    root.quit()
    root.destroy()

# Tworzenie okna Tkinter
root = tk.Tk()
root.title("Symulacja układu dynamicznego")

# Pobranie rozdzielczości ekranu
screen_width = root.winfo_screenwidth()
screen_height = root.winfo_screenheight()

# Dynamiczne skalowanie wykresów do rozdzielczości
fig_width = screen_width * 0.5 / 100  # Skalujemy wykres na 50% szerokości ekranu
fig_height = screen_height * 0.8 / 100  # Skalujemy wykres na 80% wysokości ekranu

# Tworzenie wykresu Matplotlib, dostosowanego do rozdzielczości
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(fig_width, fig_height))  # Dwa wykresy w układzie pionowym
canvas = FigureCanvasTkAgg(fig, master=root)
canvas.get_tk_widget().grid(row=0, column=0, rowspan=6)  # Ustawiamy wykres w siatce

# Tworzenie ramki na suwaki i przyciski
frame = tk.Frame(root)
frame.grid(row=0, column=1, sticky="n")  # Umieszczamy ramkę po prawej stronie wykresów

# Ustawienie suwaków jeden pod drugim w kolumnie, ale przesuwają się poziomo
slider_F = tk.Scale(frame, from_=0.1, to=10.0, resolution=0.1, orient=tk.HORIZONTAL, label="Siła F")
slider_F.grid(row=0, column=0, pady=10)

slider_M = tk.Scale(frame, from_=0.1, to=10.0, resolution=0.1, orient=tk.HORIZONTAL, label="Masa M")
slider_M.grid(row=1, column=0, pady=10)

slider_b = tk.Scale(frame, from_=0.0, to=5.0, resolution=0.1, orient=tk.HORIZONTAL, label="Tłumienie b")
slider_b.grid(row=2, column=0, pady=10)

slider_k = tk.Scale(frame, from_=0.1, to=20.0, resolution=0.1, orient=tk.HORIZONTAL, label="Sprężystość k")
slider_k.grid(row=3, column=0, pady=10)

slider_duration = tk.Scale(frame, from_=0.1, to=10.0, resolution=0.1, orient=tk.HORIZONTAL, label="Czas trwania")
slider_duration.grid(row=4, column=0, pady=10)

slider_frequency = tk.Scale(frame, from_=0.1, to=10.0, resolution=0.1, orient=tk.HORIZONTAL, label="Częstotliwość")
slider_frequency.grid(row=5, column=0, pady=10)

# Opcje wyboru rodzaju sygnału obok suwaków
radio_frame = tk.Frame(frame)  # Tworzymy nową ramkę wewnątrz głównej ramki
radio_frame.grid(row=0, column=1, rowspan=6, padx=20)  # Przyciski po prawej stronie suwaków

signal_type_var = tk.StringVar(value="prostokątny")

radio_rectangular = tk.Radiobutton(radio_frame, text="Prostokątny", variable=signal_type_var, value="prostokątny")
radio_rectangular.grid(row=0, column=0, sticky="w")

radio_triangular = tk.Radiobutton(radio_frame, text="Trójkątny", variable=signal_type_var, value="trójkątny")
radio_triangular.grid(row=1, column=0, sticky="w")

radio_harmonic = tk.Radiobutton(radio_frame, text="Harmoniczny", variable=signal_type_var, value="harmoniczny")
radio_harmonic.grid(row=2, column=0, sticky="w")

# Funkcja do obsługi suwaków i sygnałów
def on_change(event=None):
    F_const = slider_F.get()
    M = slider_M.get()
    b = slider_b.get()
    k = slider_k.get()
    duration = slider_duration.get()
    frequency = slider_frequency.get()
    signal_type = signal_type_var.get()
    update_plot(F_const, M, b, k, signal_type, duration, frequency)

# Zmieniamy zdarzenie na `"<ButtonRelease-1>"` zamiast `"<Motion>"`
slider_F.bind("<ButtonRelease-1>", on_change)
slider_M.bind("<ButtonRelease-1>", on_change)
slider_b.bind("<ButtonRelease-1>", on_change)
slider_k.bind("<ButtonRelease-1>", on_change)
slider_duration.bind("<ButtonRelease-1>", on_change)
slider_frequency.bind("<ButtonRelease-1>", on_change)
radio_rectangular.config(command=on_change)
radio_triangular.config(command=on_change)
radio_harmonic.config(command=on_change)

# Uruchomienie pętli głównej Tkinter
root.mainloop()
