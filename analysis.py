import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.optimize import curve_fit
from main import run_gibbs_ensemble
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
## Here we do the data analisys

densities,avg_v_fraction, avg_n_fraction,T = run_gibbs_ensemble()
info_equilibrio = [[],[],[],[]] #Aqui se guarda la informacion de equilibrio, cuatro entradas: densidad de la fase ,gas,densidad critica (en caso de haber),liquida y temperatura

#Primero graficamos para encontrar a que numero de ciclos se encuentra en equibilirbio
dens1, dens2 = zip(*densities)
cycles = list(range(len(densities)))

plt.plot(cycles, dens1, label="Caja 1")
plt.plot(cycles, dens2, label="Caja 2")

# Estética
plt.xlabel("Ciclo")
plt.ylabel("Densidad")
plt.title("Evolución del promedio de densidades")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("densidades_promedio.pdf", dpi=300)
plt.show()

#Usando esto encontramos el ciclo en el que se obtiene equiblibrio
ciclo_equilibrio = 10 #el 10 solo es ejemplo

#Hacemos un histograma para encontrar las densidades de equilibrio
n_gauss = 2 #es el numero de gaussianas cada pico corresponde a la densidad liquida y a la desnidad gas, pero ojo si estamos cercas de punto critico tendremos tres picos.
############### EStA FUNCIÓN ES LA QUE HACE EL AJUSTE A LAS GAUSSIANAS #############################
def ajustar_gaussianas(densidades,ciclo, n_gauss=2, bins=30, filename="ajuste_gaussianas.pdf"):
    """
    Ajusta una suma de n_gauss gaussianas al histograma de densidades.

    Args:
        densidades (list of tuple): Lista de tuplas con densidades, por ejemplo [(ρ1, ρ2), ...]
        n_gauss (int): Número de gaussianas a ajustar (2 o 3).
        bins (int): Número de bins para el histograma.
        filename (str): Nombre del archivo donde se guarda la figura.

    Returns:
        list of dicts: Parámetros de las gaussianas ajustadas (amplitud, media, sigma).
    """
    # Concatenar todas las densidades en una lista plana
    densities = densidades[ciclo:]
    all_densities = [val for pair in densities for val in pair]

    # Histograma
    counts, bin_edges = np.histogram(all_densities, bins=bins)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

    # Construir función con n gaussianas
    def multi_gaussian(x, *params):
        result = 0
        for i in range(n_gauss):
            A = params[i*3]
            mu = params[i*3 + 1]
            sigma = params[i*3 + 2]
            result += A * np.exp(-((x - mu)**2) / (2 * sigma**2))
        return result

    # Estimación inicial automática
    A0 = max(counts)
    p0 = []
    interval = (max(all_densities) - min(all_densities)) / (n_gauss + 1)
    for i in range(n_gauss):
        p0 += [A0, min(all_densities) + (i+1)*interval, 0.01]

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
        sigma = params_opt[i*3 + 2]
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
    plt.close()

    return gaussian_params
#Hacemos el ajuste
params = ajustar_gaussianas(densities,ciclo_equilibrio ,n_gauss)

for i, g in enumerate(params):
    print(f"Gaussiana {i+1}: mu = {g['mu']:.4f}, sigma = {g['sigma']:.4f}, A = {g['A']:.2f}")
#Guardamos la informacion
info_equilibrio.append((parametros,T))
