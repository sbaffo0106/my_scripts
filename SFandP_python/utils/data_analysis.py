import numpy as np 
from sklearn.utils.extmath import randomized_svd 
import matplotlib.pyplot as plt     


def tls_sf(x1, y1):
    x_centered, y_centered = x1 - np.mean(x1), y1 - np.mean(y1)
    A = np.column_stack((x_centered, y_centered))
    U, S, Vt = randomized_svd(A, n_components=2, random_state=42)
    a, b = Vt[-1]
    slope = -a / b 
    intercept = np.mean(y1) - slope * np.mean(x1)
    return slope, intercept
    
def tls_passive(x2, y2):
    x_centered2, y_centered2 = x2 - np.mean(x2), y2 - np.mean(y2)
    A2 = np.column_stack((x_centered2, y_centered2))
    U2, S2, Vt2 = randomized_svd(A2, n_components=2, random_state=42)
    a2, b2 = Vt2[-1]
    slope2 = -a2 / b2 
    intercept2 = np.mean(y2) - slope2 * np.mean(x2)
    return slope2, intercept2
    
def projection(x0, y0, slope_sf, intercept_sf, slope_passive, intercept_passive):
    log_mass_intersection = (intercept_passive - intercept_sf) / (slope_sf - slope_passive)
    SFR_intersection = slope_sf * log_mass_intersection + intercept_sf
    SFR_point = slope_sf * x0 + intercept_sf
    distanza = np.sqrt((log_mass_intersection - x0)**2 + (SFR_intersection - SFR_point)**2)
    A = 1 + slope_passive**2
    B = -2 * log_mass_intersection + 2 * intercept_passive * slope_passive - 2 * slope_passive * SFR_intersection
    C = SFR_intersection**2 - 2 * SFR_intersection * intercept_passive - distanza**2 + intercept_passive**2 + log_mass_intersection**2
    x_proj = (-B - np.sqrt(B**2 - 4 * A * C)) / (2 * A)
    y_proj = slope_passive * x_proj + intercept_passive
    return x_proj, y_proj
    
def plot_iteration(x, y, sf_filter, pass_filter, sf_tls, pass_tls, green_line, iteration):
    plt.figure(figsize=(7, 8))
    plt.scatter(x[sf_filter], y[sf_filter], c='cyan', s=2, label='Star-forming', alpha=0.5)
    plt.scatter(x[pass_filter], y[pass_filter], c='pink', s=2, label='Passive', alpha=0.5)
    contour = plt.contour(X, Y, Z, levels=6, cmap='gray')
    sf_slope, sf_intercept = sf_tls
    pass_slope, pass_intercept = pass_tls
    plt.plot(x_values2, sf_slope * x_values2 + sf_intercept, 'blue', label='TLS Star-forming')
    plt.plot(x_values2, pass_slope * x_values2 + pass_intercept, 'red', label='TLS Passive')
    green_x, green_y = green_line
    plt.plot(green_x, green_y, 'green', linewidth=2, label=f'Green Line Iter {iteration}')  
    ##print(slope_green, intercept_green)
    plt.xlim(7.5, 12)
    plt.ylim(-4.5, 2)
    plt.xlabel('M [M$_{\odot}$ h$_{70}^{-2}$]')
    plt.ylabel(r'H$\alpha$ SFR [M$_{\odot}$ yr$^{-1}$]')
    plt.legend()
    plt.title(f'Iteration {iteration}')
    plt.grid(True)
    ##plt.show()