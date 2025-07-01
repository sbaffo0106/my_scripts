import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from sklearn.utils.extmath import randomized_svd
from scipy.stats import linregress

from utils.data_processing import load_and_process_fits
from utils.data_analysis import tls_sf, tls_passive, projection
from utils.plotting import plot_results


cataid, log_mass, log_mass_err, log_sfr_halpha = load_and_process_fits('/Users/antonio/Desktop/my_scripts/data/allgalaxies.fits')

#SF mask
mask = (log_sfr_halpha - log_mass) >= -10.8
mass_filtered = log_mass[mask]
mass_err_filtered = log_mass_err[mask]
sfr_halpha_filtered = log_sfr_halpha[mask]

#mass and sfr are centered by subtracting the mean: to simplify calculations in the SVD
mass_mean = np.mean(mass_filtered)
mass_centered = mass_filtered - mass_mean
sfr_halpha_mean = np.mean(sfr_halpha_filtered)
sfr_halpha_centered = sfr_halpha_filtered - sfr_halpha_mean

#TLS is solved by performing a Singular Value Decomposition (SVD): A = USVt, where U is an
#orthogonal matrix describing the directions in the row-space of A, S a diagonal matrix containing 
#the singular values representing the variance along different directions, and Vt an 
#orthogonal matrix describing the directions in the column-space of A.
A = np.column_stack((mass_centered, sfr_halpha_centered))
U, S, Vt = randomized_svd(A, n_components=2, random_state=42) #fast approximate SVD (for large datasets)
a, b = Vt[-1]  #last row corresponds to the smallest singular value, i.e. the direction of TLS 
                     #(direction of minimal variance, which is the line of best fit in TLS)
                     #this vector is orthogonal to the TLS line in the centered space
slope_tls = -a / b #based on the relationship between the orthogonal direction vector and the TLS line
intercept_tls = sfr_halpha_mean - slope_tls * mass_mean #back to the original coordinate system
print(slope_tls,intercept_tls)

#PASSIVE mask
mask2 = (log_sfr_halpha - log_mass) < -10.8
mass_filtered2 = log_mass[mask2]
mass_err_filtered2 = log_mass_err[mask2]
sfr_halpha_filtered2 = log_sfr_halpha[mask2]
mass_mean2 = np.mean(mass_filtered2)
mass_centered2 = mass_filtered2 - mass_mean2
sfr_halpha_mean2 = np.mean(sfr_halpha_filtered2)
sfr_halpha_centered2 = sfr_halpha_filtered2 - sfr_halpha_mean2
A2 = np.column_stack((mass_centered2, sfr_halpha_centered2))
U2, S2, Vt2 = randomized_svd(A2, n_components=2, random_state=42)
a2, b2 = Vt2[-1]
slope_tls2 = -a2 / b2 
intercept_tls2 = sfr_halpha_mean2 - slope_tls2 * mass_mean2
print(slope_tls2,intercept_tls2)

#iterative method
threshold1 = 1
threshold2 = -11
max_iterations = 20
prev_green_line = None
for iteration in range(max_iterations):
    sf_filter = log_sfr_halpha >= threshold1 * log_mass + threshold2
    pass_filter = log_sfr_halpha < threshold1 * log_mass + threshold2
    if np.sum(sf_filter) == 0 or np.sum(pass_filter) == 0:
        print(f"Iterazione {iteration}: uno dei filtri è vuoto. Arresto.")
        break
    sf_tls = tls_sf(log_mass[sf_filter], log_sfr_halpha[sf_filter])
    pass_tls = tls_passive(log_mass[pass_filter], log_sfr_halpha[pass_filter])
    x_values = np.linspace(9.35, 10.75, 5)
    y_values_sf = sf_tls[0] * x_values + sf_tls[1]
    blue_x = [] 
    blue_y = [] 
    red_x = []  
    red_y = []  
    for i in range(len(x_values)):
        x0, y0 = x_values[i], y_values_sf[i] 
        x_perp, y_perp = projection(x0, y0, slope_tls, intercept_tls, slope_tls2, intercept_tls2)
        blue_x.append(x0)
        blue_y.append(y0)
        red_x.append(x_perp)
        red_y.append(y_perp)
    n_points = len(blue_x)
    x_min_all = []
    y_min_all = []
    for ii in range(n_points - 1):
        x11_blue, y11_blue = blue_x[ii], blue_y[ii]
        x12_blue, y12_blue = blue_x[ii + 1], blue_y[ii + 1]
        x11_red, y11_red = red_x[ii], red_y[ii]
        x12_red, y12_red = red_x[ii + 1], red_y[ii + 1]
        slope_sf = (y12_blue - y11_blue) / (x12_blue - x11_blue)
        intercept_sf = y11_blue - slope_sf * x11_blue
        slope_passive = (y12_red - y11_red) / (x12_red - x11_red)
        intercept_passive = y11_red - slope_passive * x11_red
        slope_1s = (y11_blue - y11_red) / (x11_blue - x11_red)
        intercept_1s = y11_blue - slope_1s * x11_blue
        slope_2s = (y12_blue - y12_red) / (x12_blue - x12_red)
        intercept_2s = y12_blue - slope_2s * x12_blue
        def inside_trapezoid(x, y):
            return (
                (y <= (slope_sf * x + intercept_sf)) and
                (y >= (slope_passive * x + intercept_passive)) and
                (y >= (slope_1s * x + intercept_1s)) and
                (y <= (slope_2s * x + intercept_2s))
            )
        mask_trap = [inside_trapezoid(x, y) for x, y in zip(log_mass, log_sfr_halpha)]
        log_mass_trap = log_mass[mask_trap]
        log_sfr_halpha_trap = log_sfr_halpha[mask_trap]
        xmin, xmax = log_mass_trap.min(), log_mass_trap.max()
        ymin, ymax = log_sfr_halpha_trap.min(), log_sfr_halpha_trap.max()
        bins = 20
        H, xedges, yedges = np.histogram2d(
            log_mass_trap, log_sfr_halpha_trap, bins=bins, range=[[xmin, xmax], [ymin, ymax]]
        )
        mask_zero = H == 0
        mask_adjacent_to_zero = np.zeros_like(H, dtype=bool)
        for row in range(H.shape[0]):
            for col in range(H.shape[1]):
                if H[row, col] == 0:
                    mask_adjacent_to_zero[max(row-2, 0):min(row+3, H.shape[0]), max(col-2, 0):min(col+3, H.shape[1])] = True
        H[mask_zero | mask_adjacent_to_zero] = 100
        min_density = np.min(H)
        min_index = np.unravel_index(np.argmin(H), H.shape)
        x_min = (xedges[min_index[0]] + xedges[min_index[0] + 1]) / 2
        y_min = (yedges[min_index[1]] + yedges[min_index[1] + 1]) / 2
        x_min_all.append(x_min)
        y_min_all.append(y_min)
    slope_green, intercept_green, _, _, _ = linregress(x_min_all, y_min_all)
    print(slope_green, intercept_green)
    x_values2 = np.linspace(7, 12, 100)
    green_line = x_values2, slope_green * x_values2 + intercept_green
    ##plot_iteration(log_mass, log_sfr_halpha, sf_filter, pass_filter, sf_tls, pass_tls, green_line, iteration)
    if prev_green_line is not None:
        x_check = np.linspace(9.35, 10.75, 5)
        y_current = np.interp(x_check, green_line[0], green_line[1])
        y_previous = np.interp(x_check, prev_green_line[0], prev_green_line[1])
        diff = np.abs(y_current - y_previous)
        print(diff)
        if np.all(diff < 0.001):
            print(f'Convergence reached after {iteration} iterations.')
            ##print(f"Slope dl: {slope_green:.2f}, Intercept dl: {intercept_green:.2f}")
            print(slope_green, intercept_green)
            break        
    prev_green_line = green_line
    threshold1 = slope_green
    threshold2 = intercept_green
sf_slope, sf_intercept = sf_tls
pass_slope, pass_intercept = pass_tls
print(f"Slope sf: {sf_slope:.2f}, Intercept sf: {sf_intercept:.2f}")
print(f"Slope passive: {pass_slope:.2f}, Intercept passive: {pass_intercept:.2f}")

plot_results(log_mass, log_sfr_halpha, sf_slope, sf_intercept, pass_slope, pass_intercept, slope_green, intercept_green)
