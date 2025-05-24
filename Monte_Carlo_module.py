#Importamos los modulos de energia
from Energy_module import pair_energy, total_energy
# =============================================================================
# Monte Carlo Move Functions
# =============================================================================
def displacement_move(x,cajas, npart, beta, vmax):
     # Energía del estado actual
    vmin = 0.01 * vmax  # esto es para evitar que una caja colapse
    vmax_abs = 10 * vmax
    enlo1 = total_energy(x,0,cajas[0][0])
    enlo2 = total_energy(x,1,cajas[1][0])

    vol1 = cajas[0][0]**3
    vol2 =  cajas[1][0]**3
    vol_diff = vol2 - vol1

    # Random walk en ln(V1/V2)
    lnvn = np.log(vol1 / vol2) + (np.random.rand() - 0.5) * vmax/5
    v1n = vol1 * np.exp(lnvn) / (1 + np.exp(lnvn))
    v2n = vol1 + vol2 - v1n
    if v1n < vmin or v2n < vmin or v1n > vmax_abs or v2n > vmax_abs: #Si una caja colapsó rechazamos
        return cajas[0][0], cajas[1][0]
    box1n = v1n**(1/3)
    box2n = v2n**(1/3)
    # Reescalado de posiciones
    for i in range(npart):
        if x[1][i] == 0:
            factor = box1n / cajas[0][0]
        else:
            factor = box2n / cajas[1][0]
        x[0][i] *= factor

    # Energía en la nueva configuración
    en1n = total_energy(x,0,box1n)
    en2n = total_energy(x,1,box2n)

    #Exponente del criterio de aceptacion
    arg1 = -beta * (en1n - enlo1) + (cajas[0][1] + 1) * np.log(v1n / vol1) / beta
    arg2 = -beta * (en2n - enlo2) + (cajas[1][1] + 1) * np.log(v2n / vol2) / beta

    # Regla de aceptación
    if np.random.rand() > np.exp(arg1 + arg2):
        # Rechazado: restaurar configuración antigua
        for i in range(npart):
            if x[1][i] == 0:
                factor = cajas[0][0] / box1n
            else:
                factor = cajas[1][0]/ box2n
            x[0][i] *= factor
        return cajas[0][0], cajas[1][0]# Retorna los valores originales

    return box1n, box2n  # Aceptado: retorna nuevos valores

def volume_move():
    return None

def transfer_move():
    return None
