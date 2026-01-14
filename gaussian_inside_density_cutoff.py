#!/usr/bin/env python3
"""
Gaussian-based cutoff determination for cryo-ET template matching scores.

This script analyzes local maximum peaks from template matching score maps,
separating them into "inside" (true positives within a mask) and "outside"
(background) populations. It fits Gaussian distributions to both populations
and determines an optimal cutoff threshold based on the inside population's
statistics.

Usage:
    python gaussian_inside_density_cutoff.py \
        --score_file SCORES.mrc \
        --mask_file MASK.mrc \
        [--min_distance 16] \
        [--bins 100] \
        [--sigma_mult 2.0] \
        [--output density_cutoff.png]
"""

import argparse
import sys
import numpy as np
import mrcfile
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.ndimage import maximum_filter
from scipy.stats import norm


def load_mrc(filepath, dtype=None, binary=False):
    """
    Load MRC file with error handling.
    
    Parameters
    ----------
    filepath : str
        Path to MRC file
    dtype : numpy dtype, optional
        Convert data to this dtype
    binary : bool
        If True, convert to binary mask (values > 0)
    
    Returns
    -------
    numpy.ndarray
        Loaded volume
    """
    try:
        with mrcfile.open(filepath, permissive=True) as mrc:
            data = mrc.data.copy()
            
        if binary:
            return data > 0
        if dtype is not None:
            return data.astype(dtype)
        return data
        
    except FileNotFoundError:
        print(f"Error: File not found: {filepath}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error loading {filepath}: {e}", file=sys.stderr)
        sys.exit(1)


def find_local_maxima(scores, min_distance):
    """
    Find local maxima in 3D score map.
    
    Parameters
    ----------
    scores : numpy.ndarray
        3D array of scores
    min_distance : int
        Minimum voxel separation between peaks
    
    Returns
    -------
    numpy.ndarray
        Boolean mask of local maxima positions
    """

    footprint = np.ones((min_distance,) * 3, dtype=bool)
    local_max = maximum_filter(scores, footprint=footprint)
    peaks_mask = scores == local_max

    
    return peaks_mask


def fit_gaussian(data):
    """
    Fit Gaussian distribution to data.
    
    Parameters
    ----------
    data : numpy.ndarray
        1D array of values
    
    Returns
    -------
    tuple
        (mean, std) of fitted Gaussian
    """
    if len(data) == 0:
        return np.nan, np.nan
    return norm.fit(data)


def plot_distributions(inside_peaks, outside_peaks, cutoff, sigma_mult, 
                       bins, output_path):
    """
    Create visualization of peak score distributions and cutoff.
    
    Parameters
    ----------
    inside_peaks : numpy.ndarray
        Peak scores inside mask
    outside_peaks : numpy.ndarray
        Peak scores outside mask
    cutoff : float
        Calculated cutoff threshold
    sigma_mult : float
        Sigma multiplier used for cutoff
    bins : int
        Number of histogram bins
    output_path : str
        Path to save figure
    """
    # Fit distributions
    mu_in, sigma_in = fit_gaussian(inside_peaks)
    mu_out, sigma_out = fit_gaussian(outside_peaks)
    
    # Create figure
    fig, ax = plt.subplots(figsize=(8, 5))
    
    # Plot histograms
    ax.hist(outside_peaks, bins=bins, density=True, alpha=0.5, 
            color='C0', edgecolor='black', linewidth=0.5,
            label=f'Outside: μ={mu_out:.3f}, σ={sigma_out:.3f}')
    ax.hist(inside_peaks, bins=bins, density=True, alpha=0.5, 
            color='C1', edgecolor='black', linewidth=0.5,
            label=f'Inside: μ={mu_in:.3f}, σ={sigma_in:.3f}')
    
    # Plot fitted Gaussian PDFs
    all_peaks = np.concatenate([inside_peaks, outside_peaks])
    x_range = np.linspace(all_peaks.min(), all_peaks.max(), 1000)
    ax.plot(x_range, norm.pdf(x_range, mu_out, sigma_out), 
            'C0--', linewidth=2, alpha=0.8, label='Outside fit')
    ax.plot(x_range, norm.pdf(x_range, mu_in, sigma_in), 
            'C1--', linewidth=2, alpha=0.8, label='Inside fit')
    
    # Mark cutoff threshold
    ymax = ax.get_ylim()[1]
    ax.axvline(cutoff, ymin=0, ymax=0.05, color='red', linewidth=4)
    ax.plot([], [], color='red', linewidth=4,
            label=f'Cutoff: μ−{sigma_mult:.1f}σ = {cutoff:.3f}')
    
    # Formatting
    ax.set_xlabel('Peak Score (CC)', fontsize=12)
    ax.set_ylabel('Probability Density', fontsize=12)
    ax.set_title(f'Score Distribution Analysis (μ−{sigma_mult:.1f}σ cutoff)', 
                 fontsize=13, fontweight='bold')
    ax.legend(loc='best', fontsize=9, framealpha=0.9)
    ax.grid(True, alpha=0.3, linestyle='--')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Plot {output_path} saved")


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--score_file", required=True,
        help="MRC file containing template matching scores"
    )
    parser.add_argument(
        "--mask_file", required=True,
        help="MRC file containing binary mask (0=background, 1=target region)"
    )
    parser.add_argument(
        "--min_distance", type=int, default=16,
        help="Minimum voxel separation for local maxima (default: 16)"
    )
    parser.add_argument(
        "--bins", type=int, default=100,
        help="Number of histogram bins (default: 100)"
    )
    parser.add_argument(
        "--sigma_mult", type=float, default=2.0,
        help="Sigma multiplier for cutoff threshold (default: 2.0)"
    )
    parser.add_argument(
        "--output", default="density_cutoff.png",
        help="Output PNG filename (default: density_cutoff.png)"
    )
    
    args = parser.parse_args()
    
    # Load data
    print(f"Loading score map: {args.score_file}")
    scores = load_mrc(args.score_file, dtype=np.float32)
    
    print(f"Loading mask: {args.mask_file}")
    mask = load_mrc(args.mask_file, binary=True)
    
    # Validate dimensions
    if scores.shape != mask.shape:
        print(f"Error: Score map {scores.shape} and mask {mask.shape} "
              f"dimensions do not match", file=sys.stderr)
        sys.exit(1)
    
    # Find local maxima
    print(f"Finding local maxima (min_distance={args.min_distance})...")
    peaks_mask = find_local_maxima(scores, args.min_distance)
    
    # Extract peak values
    all_peaks = scores[peaks_mask]
    inside_peaks = all_peaks[mask[peaks_mask]]
    outside_peaks = all_peaks[~mask[peaks_mask]]
    

    if len(inside_peaks) == 0:
        print("Error: No peaks found inside mask!", file=sys.stderr)
        sys.exit(1)
    
    # Calculate statistics
    mu_in, sigma_in = fit_gaussian(inside_peaks)
    mu_out, sigma_out = fit_gaussian(outside_peaks)
    
    print(f"\nDistribution statistics:")
    print(f"  Inside:  μ = {mu_in:.4f}, σ = {sigma_in:.4f}")
    print(f"  Outside: μ = {mu_out:.4f}, σ = {sigma_out:.4f}")
    
    # Calculate cutoff
    cutoff = mu_in - args.sigma_mult * sigma_in
    
    print(f"\n{'='*50}")
    print(f"RECOMMENDED CUTOFF: {cutoff:.4f}")
    print(f"{'='*50}\n")
    
    # Create visualization
    print("Generating plot...")
    plot_distributions(inside_peaks, outside_peaks, cutoff, 
                      args.sigma_mult, args.bins, args.output)
    
    print("\nDone!")


if __name__ == '__main__':
    main()