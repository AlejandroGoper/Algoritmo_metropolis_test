import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.optimize import curve_fit
from main import run_gibbs_ensemble
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
## Here we do the data analisys


#Usando esto encontramos el ciclo en el que se obtiene equiblibrio
ciclo_equilibrio = 10 #el 10 solo es ejemplo

#Hacemos un histograma para encontrar las densidades de equilibrio
n_gauss = 2 #es el numero de gaussianas cada pico corresponde a la densidad liquida y a la desnidad gas, pero ojo si estamos cercas de punto critico tendremos tres picos.
############### EStA FUNCIÓN ES LA QUE HACE EL AJUSTE A LAS GAUSSIANAS #############################
def ajustar_gaussianas(densidades, ciclo, n_gauss=2, bins=50, filename="ajuste_gaussianas.pdf",dist = 3):
    """
    Ajusta una suma de n_gauss gaussianas al histograma de densidades.

    Args:
        densidades (list of tuple): Lista de tuplas con densidades, por ejemplo [(ρ1, ρ2), ...]
        ciclo (int): Índice a partir del cual se usan las densidades.
        n_gauss (int): Número de gaussianas a ajustar (2 o 3).
        bins (int): Número de bins para el histograma.
        filename (str): Nombre del archivo donde se guarda la figura.

    Returns:
        list of dicts: Parámetros de las gaussianas ajustadas (amplitud, media, sigma).
    """
    import math

    # Truncar densidades y aplanar
    densities = densidades[ciclo:]
    all_densities = [val for pair in densities for val in pair]

    # Histograma
    counts, bin_edges = np.histogram(all_densities, bins=bins)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

    # Encontrar picos del histograma
    peaks, _ = find_peaks(counts, distance=dist)  # Ajusta 'distance' si hay picos muy cercanos
        
    # Ordenar los picos por altura descendente y tomar los n_gauss más prominentes
    sorted_peaks = sorted(peaks, key=lambda x: counts[x], reverse=True)[:n_gauss]
    sorted_peaks = sorted(sorted_peaks)  # ordenarlos por posición para consistencia

    # Estimación inicial basada en los picos
    p0 = []
    for peak in sorted_peaks:
        A = counts[peak]
        mu = bin_centers[peak]
        sigma = 0.01  # o estima localmente usando los bins si quieres más precisión
        p0 += [A, mu, sigma]

    # Si no hay suficientes picos, completa con valores distribuidos
    while len(p0) < 3 * n_gauss:
        p0 += [max(counts), np.mean(all_densities), 0.01]

    # Construir función con n gaussianas
    def multi_gaussian(x, *params):
        result = 0
        for i in range(n_gauss):
            A = params[i*3]
            mu = params[i*3 + 1]
            sigma = params[i*3 + 2]
            result += A * np.exp(-((x - mu)**2) / (2 * sigma**2))
        return result

    # Ajuste
    try:
        params_opt, _ = curve_fit(multi_gaussian, bin_centers, counts, p0=p0)
    except RuntimeError as e:
        print("Error en el ajuste:", e)
        return []

    # Extraer parámetros
    gaussian_params = []
    for i in range(n_gauss):
        A = params_opt[i*3]
        mu = params_opt[i*3 + 1]
        sigma = abs(params_opt[i*3 + 2])  # asegurar positivo
        gaussian_params.append({'A': A, 'mu': mu, 'sigma': sigma})

    # Evaluar curva ajustada
    x_fit = np.linspace(min(bin_centers), max(bin_centers), 1000)
    y_fit = multi_gaussian(x_fit, *params_opt)

    # Graficar
    plt.hist(all_densities, bins=bins, alpha=0.5, edgecolor='black', label="Datos")
    plt.plot(x_fit, y_fit, 'r-', label=f"Ajuste ({n_gauss} gaussianas)", linewidth=2)
    plt.xlabel("Densidad")
    plt.ylabel("Frecuencia")
    plt.title(f"Ajuste de {n_gauss} gaussianas al histograma")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(filename, dpi=300)
    plt.show()
    plt.close()

    return gaussian_params


