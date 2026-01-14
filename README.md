# Cryo-ET Template Matching Analysis Tools

A suite of Python utilities for analyzing template matching results and particle distributions in cryo-electron tomography (cryo-ET) data.

## Overview

This toolkit provides two complementary analysis workflows:

1. **Cutoff Determination** (`gaussian_inside_density_cutoff.py`): Statistically determine optimal cross-correlation (CC) score thresholds for particle picking using a binary mask
2. **Spatial Analysis** (`nearest_neighbor.py`): Analyze spatial distribution patterns of picked particles

## Tools

### 1. Gaussian-Based Cutoff Determination

Analyzes template matching score maps to determine an optimal threshold for particle selection using Gaussian distribution fitting.

#### Features

- Identifies local maxima in 3D score maps
- Separates peaks into "inside" (true positives) and "outside" (background) populations using a tomogram mask
- Fits Gaussian distributions to both populations
- Calculates cutoff threshold: **μ_inside - N×σ_inside**

#### Usage

```bash
python gaussian_inside_density_cutoff.py \
    --score_file template_matching_scores.mrc \
    --mask_file reference_mask.mrc \
    --min_distance 28 \
    --sigma_mult 2.0 \
    --bins 100 \
    --output cutoff_analysis.png
```

#### Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--score_file` | (required) | MRC file containing template matching scores |
| `--mask_file` | (required) | Binary MRC mask **(1=target region, 0=background)** |
| `--min_distance` | 23 | Minimum **voxel** separation between local maxima |
| `--sigma_mult` | 2.0 | Sigma multiplier for cutoff |
| `--bins` | 100 | Number of histogram bins for visualization |
| `--output` | `density_cutoff.png` | Output figure filename |

#### Output

The script outputs:
- **Console**: Statistical summary and recommended cutoff value
- **PNG Figure**: Density histograms with fitted Gaussians and marked cutoff threshold


### 2. Nearest Neighbor Spatial Analysis

Analyzes the spatial distribution of particles using **RELION5 STAR** files.

#### Features

- **Nearest neighbor distance distributions**
- **Kernel density estimation**
- **Neighbor counts within spherical distances**
- **Radial distribution function (RDF)**


**Note**: Edit the following parameters directly in the script:

```python
star_path   = "path/to/particles.star"   # Path to your STAR file
k_max       = 8                          # Number of nearest neighbors
xlim        = (100, 300)                 # X-axis limits for distance plots (Å)
rdf_r_max   = 800                        # Max radius for RDF (Å)
rdf_dr      = 20                         # RDF bin width (Å)
radii_nm    = [15, 20, 25, 30, 35, 40]   # Radii for neighbor counts (nm)
```

#### Required STAR File Columns

The script expects the following columns in your STAR file:
- `rlnTomoName`: Tomogram identifier
- `rlnCenteredCoordinateXAngst`: X coordinate (Ångströms)
- `rlnCenteredCoordinateYAngst`: Y coordinate (Ångströms)
- `rlnCenteredCoordinateZAngst`: Z coordinate (Ångströms)

#### Output

For each tomogram, generates a 4-panel figure:

1. **Histogram**: Distance distributions nearest neighbors
2. **Density Plot**: KDE-smoothed distributions with peak positions
3. **Neighbor Counts**: Boxplots showing particle counts within specified radii
4. **RDF**: Radial distribution function revealing spatial patterns


### Requirements

```bash
pip install numpy scipy matplotlib mrcfile starfile
```

## License

MIT