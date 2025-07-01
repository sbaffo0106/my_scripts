import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from matplotlib.colors import LinearSegmentedColormap

def plot_results(log_mass, log_sfr_halpha, sf_slope, sf_intercept, pass_slope, pass_intercept, slope_green, intercept_green):
    colors = [(0.5, 0.5, 0.5), (0, 0, 0)]
    custom_cmap = LinearSegmentedColormap.from_list("custom_gray", colors)
    values = np.vstack([log_mass, log_sfr_halpha])
    kde = gaussian_kde(values, bw_method=0.06)
    xmin, xmax = log_mass.min(), log_mass.max()
    ymin, ymax = log_sfr_halpha.min(), log_sfr_halpha.max()
    xi, yi = np.linspace(xmin, xmax, 100), np.linspace(ymin, ymax, 100)
    X, Y = np.meshgrid(xi, yi)
    positions = np.vstack([X.ravel(), Y.ravel()])
    Z = kde(positions).reshape(X.shape)

    plt.figure(figsize=(7, 7))
    plt.contour(X, Y, Z, levels=8, cmap=custom_cmap)
    x_range = np.linspace(7.5, 12, 100)
    y_sf = sf_slope * x_range + sf_intercept
    plt.plot(x_range, y_sf, '--', color='blue', linewidth=2.5, label='SF fit')
    y_passive = pass_slope * x_range + pass_intercept
    plt.plot(x_range, y_passive, '--', color='red', linewidth=2.5, label='Passive fit')
    plt.scatter(log_mass, log_sfr_halpha, s=2, color='gray')
    plt.plot(x_range, slope_green * x_range + intercept_green, '-', color='green', linewidth=2.5, label='demarcation line')
    above_green = log_sfr_halpha > (slope_green * log_mass + intercept_green)
    below_green = ~above_green
    plt.scatter(log_mass[above_green], log_sfr_halpha[above_green], s=2, color='cyan')
    plt.scatter(log_mass[below_green], log_sfr_halpha[below_green], s=2, color='pink')
    plt.xlim(7.5, 12)
    plt.ylim(-4, 2)
    plt.xlabel(r'$\log\, \mathrm{M}\ [\mathrm{M}_\odot\, h_{70}^{-2}]$', fontsize=18)
    plt.ylabel(r'$\log \mathrm{SFR}_{\mathrm{H}\alpha}$ [M$_\odot$ yr$^{-1}$]', fontsize=18)
    plt.xticks([8, 9, 10, 11, 12], fontsize=14)
    plt.yticks(fontsize=14)
    plt.legend(fontsize=13)
    plt.show()