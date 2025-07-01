#Iterative Classification of Star-Forming and Passive Galaxies
This project implements an iterative method to separate star-forming and passive galaxies based on their position in the stellar mass–SFR plane, using Hα-derived star formation rates and stellar masses.

##Background
Traditional methods to distinguish between blue/star-forming and red/passive galaxies often rely on arbitrary cuts in color–magnitude or SFR–mass diagrams. However, at low stellar masses, dust attenuation and diverse star formation histories complicate this separation. More sophisticated probabilistic models (e.g., Taylor et al. 2015) offer better accuracy but are computationally intensive and tailored for specific datasets.

An alternative method was proposed by Davies et al. (2019), who defined the boundary between star-forming and passive populations as the ridge line of minimum galaxy density in the SFR–mass plane, determined from perpendicular slices to two independent regression fits.

##Method
This project follows a similar approach to Davies et al. (2019) but removes the need for an initial arbitrary separation. Instead, it uses an **iterative process**:
1. Perform two Total Least Squares (TLS) fits on preliminary star-forming and passive subsets.
2. Define a dividing line based on the local minimum density between these populations.
3. Re-classify galaxies based on their position relative to this new boundary.
4. Repeat until the dividing line converges (tolerance: 0.01 dex).

The result is a self-consistent and data-driven separation between the two populations.

##Project structure
- `main.py`: Runs the full analysis pipeline.
- `utils/data_processing.py`: Loads and processes input FITS data.
- `utils/data_analysis.py`: Performs TLS fitting, projections, and iteration logic.
- `utils/plotting.py`: Generates plots of the classification and density contours.

## Requirements
Install the dependencies with:

```bash
pip install -r requirements.txt

##License
This project is intended for academic and educational use. For full terms and conditions, please refer to the LICENSE file.

##Author
Antonio Sbaffoni – 2025
