#!/usr/bin/env python3

import math
from typing import List, Tuple, Dict, Any
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Constants
mepsilon = 1e-10

def MapIndex(ir1: int, ir2: int) -> int:
    """
    Map 2D indices to 1D index using symmetry.
    This function maps (ir1, ir2) to a 1D index, assuming ir1 >= ir2.
    """
    return ir1 * (ir1 + 1) // 2 + ir2

def CosThetaMin(rmax: float, r1: float, r2: float) -> float:
    """
    Calculate minimum cosine of theta.
    """
    return max(-1.0, (r1**2 + r2**2 - rmax**2) / (2 * r1 * r2))

def CosThetaMax(rmin: float, r1: float, r2: float) -> float:
    """
    Calculate maximum cosine of theta.
    """
    return min(1.0, (r1**2 + r2**2 - rmin**2) / (2 * r1 * r2))

def ReadTBF(fname: str) -> Dict[str, Any]:
    """
    Read TBF (Three-Body Function) data from file.
    
    Args:
        fname: Path to the input file
        
    Returns:
        Dictionary containing the TBF data structure
    """
    # Read file and parse initial parameters
    with open(fname, 'r') as fid:
        lines = fid.readlines()
    
    # Parse header information
    nsample = int(lines[0].strip())
    ntype = int(lines[1].strip())
    rmin, rmax, dr = map(float, lines[2].strip().split())
    
    # Initialize TBF structure
    tbf = {
        'nsample': nsample,
        'ntype': ntype,
        'rmin': rmin,
        'rmax': rmax,
        'dr': dr
    }
    
    # Calculate number of radial bins
    nr = math.ceil((rmax - rmin) / dr)
    
    # Calculate number of bins using symmetry (r1 >= r2)
    nbin = nr * (nr + 1) // 2
    
    # Initialize arrays
    tmins = [0.0] * nbin
    tmaxs = [0.0] * nbin
    dt = [0.0] * nbin
    nt = [0] * nbin
    
    # Initialize 3D function array using nested lists
    # f3s[i][j][k] where i is type index, j is bin index, k is theta index
    f3s = []
    for i in range(ntype * ntype * ntype):
        f3s.append([None] * nbin)
    
    # Parse the data section (starting from line 4)
    data_line_idx = 3
    data_values = []
    
    # Collect all data values
    for i in range(data_line_idx, len(lines)):
        line_values = list(map(float, lines[i].strip().split()))
        data_values.extend(line_values)
    
    # Process each bin and type combination
    data_idx = 0
    for ir1 in range(nr):
        for ir2 in range(ir1 + 1):  # Only consider r1 >= r2
            idx = MapIndex(ir1, ir2)
            r1 = rmin + dr * ir1
            r2 = rmin + dr * ir2
            
            # Calculate theta bounds
            tmins[idx] = CosThetaMin(rmax, r1, r2)
            tmaxs[idx] = CosThetaMax(rmin, r1, r2)
            
            # Calculate r3 bounds
            r3min = max(rmin, r1 - r2 - dr)
            r3max = min(rmax, r1 + r2 + 2 * dr)
            
            # Calculate number of theta bins
            nt[idx] = 2 * math.ceil((r3max - r3min) / dr)
            dt[idx] = (tmaxs[idx] - tmins[idx] + mepsilon) / nt[idx]
            
            # Allocate and fill f3s for each type combination
            for type_idx in range(ntype * ntype * ntype):
                f3s[type_idx][idx] = [0.0] * nt[idx]
                
                # Fill with data values
                for it in range(nt[idx]):
                    if data_idx < len(data_values):
                        f3s[type_idx][idx][it] = data_values[data_idx]
                        data_idx += 1
    
    # Complete the TBF structure
    tbf['nr'] = nr
    tbf['f3s'] = f3s
    tbf['vols'] = None  # Not implemented in original
    tbf['nbin'] = nbin
    tbf['tmins'] = tmins
    tbf['tmaxs'] = tmaxs
    tbf['nt'] = nt
    tbf['dt'] = dt
    
    return tbf

def print_tbf_info(tbf: Dict[str, Any]) -> None:
    """
    Print information about the TBF structure for debugging.
    """
    print(f"nsample: {tbf['nsample']}")
    print(f"ntype: {tbf['ntype']}")
    print(f"rmin: {tbf['rmin']}")
    print(f"rmax: {tbf['rmax']}")
    print(f"dr: {tbf['dr']}")
    print(f"nr: {tbf['nr']}")
    print(f"nbin: {tbf['nbin']}")
    print(f"Number of type combinations: {len(tbf['f3s'])}")
    if tbf['f3s']:
        print(f"Number of bins for first type: {len(tbf['f3s'][0])}")
        if tbf['f3s'][0][0]:
            print(f"Number of theta values for first bin: {len(tbf['f3s'][0][0])}")

def plot_tbf_3d(tbf: Dict[str, Any], type_idx: int = 0, save_plot: bool = False, 
                filename: str = "tbf_3d_plot.png", max_points: int = 50000, 
                plot_type: str = 'scatter', threshold: float = None) -> None:
    """
    Create a 3D plot of the TBF data using r1, r2, and r3 as axes.
    Optimized for large datasets with multiple visualization options.
    
    Args:
        tbf: The TBF data structure returned by ReadTBF
        type_idx: Index of the type combination to plot (default: 0)
        save_plot: Whether to save the plot to a file (default: False)
        filename: Filename to save the plot (default: "tbf_3d_plot.png")
        max_points: Maximum number of points to plot (default: 50000)
        plot_type: Type of plot ('scatter', 'surface', 'contour') (default: 'scatter')
        threshold: Minimum absolute value to include in plot. If None, auto-calculated as 1% of max value.
    """
    # Extract parameters
    rmin = tbf['rmin']
    rmax = tbf['rmax']
    dr = tbf['dr']
    nr = tbf['nr']
    nbin = tbf['nbin']
    f3s = tbf['f3s'][type_idx]
    tmins = tbf['tmins']
    tmaxs = tbf['tmaxs']
    nt = tbf['nt']
    dt = tbf['dt']
    
    # Prepare data for 3D plotting
    r1_values = []
    r2_values = []
    r3_values = []
    f3_values = []
    
    # Collect all values first to calculate threshold if needed
    all_values = []
    point_data = []  # Store (r1, r2, r3, f3) tuples
    
    # Iterate through all bins to collect data
    for ir1 in range(nr):
        for ir2 in range(ir1 + 1):  # Only consider r1 >= r2
            idx = MapIndex(ir1, ir2)
            r1 = rmin + dr * ir1
            r2 = rmin + dr * ir2
            
            # Calculate r3 bounds for this (r1, r2) pair
            r3min = max(rmin, r1 - r2 - dr)
            r3max = min(rmax, r1 + r2 + 2 * dr)
            
            # Calculate theta values and corresponding r3 values
            if nt[idx] > 0:
                for it in range(nt[idx]):
                    theta = tmins[idx] + it * dt[idx]
                    # Convert cosine to angle and calculate r3
                    cos_theta = theta
                    r3_squared = r1**2 + r2**2 - 2 * r1 * r2 * cos_theta
                    if r3_squared >= 0:
                        r3 = math.sqrt(r3_squared)
                        # Only include points within valid r3 range
                        if r3min <= r3 <= r3max:
                            f3_value = f3s[idx][it]
                            all_values.append(abs(f3_value))
                            point_data.append((r1, r2, r3, f3_value))
    
    # Determine threshold
    if threshold is None:
        if all_values:
            # Auto-calculate threshold as 1% of maximum absolute value
            max_abs_value = max(all_values)
            threshold = 0.01 * max_abs_value
        else:
            threshold = 0.0
    
    # Filter points based on threshold
    for r1, r2, r3, f3_value in point_data:
        if abs(f3_value) >= threshold:
            r1_values.append(r1)
            r2_values.append(r2)
            r3_values.append(r3)
            f3_values.append(f3_value)
    
    print(f"Threshold applied: {threshold:.6f}")
    print(f"Total points before filtering: {len(point_data)}")
    print(f"Points after threshold filtering: {len(r1_values)}")
    
    # Create 3D plot
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')
    
    if plot_type == 'scatter':
        # Optimized scatter plot for large datasets
        if len(r1_values) > max_points:
            # Randomly sample points to improve performance
            indices = np.random.choice(len(r1_values), max_points, replace=False)
            r1_plot = [r1_values[i] for i in indices]
            r2_plot = [r2_values[i] for i in indices]
            r3_plot = [r3_values[i] for i in indices]
            f3_plot = [f3_values[i] for i in indices]
            point_size = 10  # Smaller points for better performance
        else:
            r1_plot, r2_plot, r3_plot, f3_plot = r1_values, r2_values, r3_values, f3_values
            point_size = 20
        
        # Create scatter plot with optimized settings
        scatter = ax.scatter(r1_plot, r2_plot, r3_plot, c=f3_plot, cmap='YlGn', 
                           alpha=0.6, s=point_size, marker='o', linewidth=0)
        
    elif plot_type == 'surface':
        # Create surface plot using grid data
        # This is more efficient for dense data
        r1_grid = np.linspace(rmin, rmax, int(np.sqrt(len(r1_values))))
        r2_grid = np.linspace(rmin, rmax, int(np.sqrt(len(r2_values))))
        R1, R2 = np.meshgrid(r1_grid, r2_grid)
        
        # Interpolate data onto the grid
        from scipy.interpolate import griddata
        points = np.column_stack((r1_values, r2_values, r3_values))
        values = f3_values
        
        # Create a grid for r3 values
        R3 = np.zeros_like(R1)
        F3 = np.zeros_like(R1)
        
        for i in range(R1.shape[0]):
            for j in range(R1.shape[1]):
                # Find nearest points for interpolation
                r1_target = R1[i,j]
                r2_target = R2[i,j]
                
                # Simple nearest neighbor interpolation
                distances = np.sqrt((np.array(r1_values) - r1_target)**2 + 
                                  (np.array(r2_values) - r2_target)**2)
                if len(distances) > 0:
                    min_idx = np.argmin(distances)
                    R3[i,j] = r3_values[min_idx]
                    F3[i,j] = f3_values[min_idx]
        
        # Create surface plot
        surf = ax.plot_surface(R1, R2, R3, facecolors=plt.cm.YlGn(F3/np.max(F3)), 
                             alpha=0.8, linewidth=0, antialiased=True)
        
    elif plot_type == 'contour':
        # Create contour plot
        # This is most efficient for very large datasets
        r1_unique = sorted(list(set(r1_values)))
        r2_unique = sorted(list(set(r2_values)))
        
        if len(r1_unique) > 100 or len(r2_unique) > 100:
            # Downsample for performance
            r1_unique = r1_unique[::max(1, len(r1_unique)//50)]
            r2_unique = r2_unique[::max(1, len(r2_unique)//50)]
        
        R1, R2 = np.meshgrid(r1_unique, r2_unique)
        R3 = np.zeros_like(R1)
        F3 = np.zeros_like(R1)
        
        # Interpolate values
        from scipy.interpolate import griddata
        points = np.column_stack((r1_values, r2_values))
        r3_interp = griddata(points, r3_values, (R1, R2), method='linear', fill_value=0)
        f3_interp = griddata(points, f3_values, (R1, R2), method='linear', fill_value=0)
        
        # Create contour plot
        contour = ax.contourf(R1, R2, r3_interp, levels=20, cmap='YlGn', alpha=0.8)
        ax.contour(R1, R2, r3_interp, levels=20, colors='black', alpha=0.3, linewidths=0.5)
    
    else:
        raise ValueError(f"Unknown plot_type: {plot_type}. Use 'scatter', 'surface', or 'contour'")
    
    # Customize the plot
    ax.set_xlabel('r1 (Å)', fontsize=12, labelpad=10)
    ax.set_ylabel('r2 (Å)', fontsize=12, labelpad=10)
    ax.set_zlabel('r3 (Å)', fontsize=12, labelpad=10)
    ax.set_title(f'3D Visualization\nUsing r1, r2, r3 as axes\nPlot type: {plot_type}', 
                fontsize=14, pad=20)
    
    # Add colorbar
    if plot_type == 'scatter':
        plt.colorbar(scatter, ax=ax, shrink=0.5, aspect=5, label='TBF Value')
    elif plot_type == 'surface':
        m = plt.cm.ScalarMappable(cmap=plt.cm.YlGn)
        m.set_array(f3_values)
        plt.colorbar(m, ax=ax, shrink=0.5, aspect=5, label='TBF Value')
    else:  # contour
        plt.colorbar(contour, ax=ax, shrink=0.5, aspect=5, label='TBF Value')
    
    # Set equal aspect ratio if possible
    max_range = max(rmax - rmin, rmax - rmin, rmax - rmin)
    mid_x = (rmin + rmax) / 2
    mid_y = (rmin + rmax) / 2
    mid_z = (rmin + rmax) / 2
    
    ax.set_xlim(mid_x - max_range/2, mid_x + max_range/2)
    ax.set_ylim(mid_y - max_range/2, mid_y + max_range/2)
    ax.set_zlim(mid_z - max_range/2, mid_z + max_range/2)
    
    plt.tight_layout()
    
    if save_plot:
        plt.savefig(filename, dpi=300, bbox_inches='tight', 
                   facecolor='white', edgecolor='none')
        print(f"Plot saved as {filename}")
    
    plt.show()

def plot_tbf_3d_performance(tbf: Dict[str, Any], type_idx: int = 0, save_plot: bool = False, 
                           filename: str = "tbf_3d_performance.png", threshold: float = None) -> None:
    """
    Create a performance-optimized 3D plot for very large datasets.
    Uses hexagonal binning and density visualization.
    
    Args:
        tbf: The TBF data structure returned by ReadTBF
        type_idx: Index of the type combination to plot (default: 0)
        save_plot: Whether to save the plot to a file (default: False)
        filename: Filename to save the plot (default: "tbf_3d_performance.png")
        threshold: Minimum absolute value to include in plot. If None, auto-calculated as 1% of max value.
    """
    # Extract parameters
    rmin = tbf['rmin']
    rmax = tbf['rmax']
    dr = tbf['dr']
    nr = tbf['nr']
    f3s = tbf['f3s'][type_idx]
    tmins = tbf['tmins']
    tmaxs = tbf['tmaxs']
    nt = tbf['nt']
    dt = tbf['dt']
    
    # Prepare data for 3D plotting
    r1_values = []
    r2_values = []
    r3_values = []
    f3_values = []
    
    # Collect all values first to calculate threshold if needed
    all_values = []
    point_data = []  # Store (r1, r2, r3, f3) tuples
    
    # Iterate through all bins to collect data
    for ir1 in range(nr):
        for ir2 in range(ir1 + 1):  # Only consider r1 >= r2
            idx = MapIndex(ir1, ir2)
            r1 = rmin + dr * ir1
            r2 = rmin + dr * ir2
            
            # Calculate r3 bounds for this (r1, r2) pair
            r3min = max(rmin, r1 - r2 - dr)
            r3max = min(rmax, r1 + r2 + 2 * dr)
            
            # Calculate theta values and corresponding r3 values
            if nt[idx] > 0:
                for it in range(nt[idx]):
                    theta = tmins[idx] + it * dt[idx]
                    cos_theta = theta
                    r3_squared = r1**2 + r2**2 - 2 * r1 * r2 * cos_theta
                    if r3_squared >= 0:
                        r3 = math.sqrt(r3_squared)
                        # Only include points within valid r3 range
                        if r3min <= r3 <= r3max:
                            f3_value = f3s[idx][it]
                            all_values.append(abs(f3_value))
                            point_data.append((r1, r2, r3, f3_value))
    
    # Determine threshold
    if threshold is None:
        if all_values:
            # Auto-calculate threshold as 1% of maximum absolute value
            max_abs_value = max(all_values)
            threshold = 0.01 * max_abs_value
        else:
            threshold = 0.0
    
    # Filter points based on threshold
    for r1, r2, r3, f3_value in point_data:
        if abs(f3_value) >= threshold:
            r1_values.append(r1)
            r2_values.append(r2)
            r3_values.append(r3)
            f3_values.append(f3_value)
    
    print(f"Performance plot - Threshold applied: {threshold:.6f}")
    print(f"Performance plot - Total points before filtering: {len(point_data)}")
    print(f"Performance plot - Points after threshold filtering: {len(r1_values)}")
    
    # Create performance-optimized plot
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')
    
    # Convert to numpy arrays for faster processing
    r1_array = np.array(r1_values)
    r2_array = np.array(r2_values)
    r3_array = np.array(r3_values)
    f3_array = np.array(f3_values)
    
    # Create a 3D histogram/binning approach for very large datasets
    n_bins = 50  # Adjust based on performance needs
    
    # Create bins
    r1_bins = np.linspace(rmin, rmax, n_bins)
    r2_bins = np.linspace(rmin, rmax, n_bins)
    r3_bins = np.linspace(rmin, rmax, n_bins)
    
    # Digitize the data
    r1_digitized = np.digitize(r1_array, r1_bins) - 1
    r2_digitized = np.digitize(r2_array, r2_bins) - 1
    r3_digitized = np.digitize(r3_array, r3_bins) - 1
    
    # Create a 3D histogram
    hist_3d, edges = np.histogramdd([r1_array, r2_array, r3_array], 
                                   bins=[r1_bins, r2_bins, r3_bins], 
                                   weights=f3_array)
    
    # Find non-zero bins for plotting
    non_zero_indices = np.where(hist_3d > 0)
    
    if len(non_zero_indices[0]) > 0:
        # Create coordinates for non-zero bins
        x_coords = (edges[0][non_zero_indices[0]] + edges[0][non_zero_indices[0] + 1]) / 2
        y_coords = (edges[1][non_zero_indices[1]] + edges[1][non_zero_indices[1] + 1]) / 2
        z_coords = (edges[2][non_zero_indices[2]] + edges[2][non_zero_indices[2] + 1]) / 2
        values = hist_3d[non_zero_indices]
        
        # Create scatter plot with binned data
        scatter = ax.scatter(x_coords, y_coords, z_coords, c=values, cmap='YlGn', 
                           alpha=0.8, s=50, marker='s')
        
        # Customize the plot
        ax.set_xlabel('r1 (Å)', fontsize=12, labelpad=10)
        ax.set_ylabel('r2 (Å)', fontsize=12, labelpad=10)
        ax.set_zlabel('r3 (Å)', fontsize=12, labelpad=10)
        ax.set_title(f'3D Binned Visualization\nUsing r1, r2, r3 as axes\n{len(non_zero_indices[0])} bins', 
                    fontsize=14, pad=20)
        
        plt.colorbar(scatter, ax=ax, shrink=0.5, aspect=5, label='TBF Value')
    
    else:
        ax.text2D(0.5, 0.5, "No data to display", transform=ax.transAxes, fontsize=20)
    
    plt.tight_layout()
    
    if save_plot:
        plt.savefig(filename, dpi=300, bbox_inches='tight', 
                   facecolor='white', edgecolor='none')
        print(f"Performance plot saved as {filename}")
    
    plt.show()

def plot_tbf_radius_theta(tbf: Dict[str, Any], type_idx: int = 0, save_plot: bool = False, 
                         filename: str = "tbf_radius_theta.png") -> None:
    """
    Create a 2D plot showing radius-theta correlation where r1=r2=r.
    Theta is the angle between r1 and r2 vectors.
    
    Args:
        tbf: The TBF data structure returned by ReadTBF
        type_idx: Index of the type combination to plot (default: 0)
        save_plot: Whether to save the plot to a file (default: False)
        filename: Filename to save the plot (default: "tbf_radius_theta.png")
    """
    # Extract parameters
    rmin = tbf['rmin']
    rmax = tbf['rmax']
    dr = tbf['dr']
    nr = tbf['nr']
    f3s = tbf['f3s'][type_idx]
    tmins = tbf['tmins']
    tmaxs = tbf['tmaxs']
    nt = tbf['nt']
    dt = tbf['dt']
    
    # Prepare data for radius-theta plotting
    r_values = []
    theta_values = []
    f3_values = []
    
    # Iterate through diagonal bins where r1 = r2 = r
    for ir in range(nr):
        # For r1 = r2, we use the diagonal where ir1 = ir2 = ir
        # MapIndex(ir, ir) gives us the bin index for r1 = r2 = r
        idx = MapIndex(ir, ir)
        r = rmin + dr * ir
        
        # Calculate theta bounds for this radius
        if nt[idx] > 0:
            for it in range(nt[idx]):
                theta = tmins[idx] + it * dt[idx]
                # Convert cosine to angle in degrees for better visualization
                angle_rad = math.acos(max(-1.0, min(1.0, theta)))
                angle_deg = math.degrees(angle_rad)
                
                # Only include valid angles and non-zero TBF values
                if 0 <= angle_deg <= 180 and f3s[idx][it] != 0:
                    r_values.append(r)
                    theta_values.append(angle_deg)
                    f3_values.append(f3s[idx][it])
    
    # Create 2D plot
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Create scatter plot with color representing TBF value
    scatter = ax.scatter(r_values, theta_values, c=f3_values, cmap='YlGn', 
                        alpha=0.7, s=30, edgecolors='none')
    
    # Customize the plot
    ax.set_xlabel('Radius r (Å)', fontsize=14, labelpad=10)
    ax.set_ylabel('Angle θ between r₁ and r₂ (degrees)', fontsize=14, labelpad=10)
    ax.set_title(f'Radius-Angle Correlation (r₁ = r₂ = r)', fontsize=16, pad=20)
    
    # Add colorbar
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label('TBF Value', fontsize=12)
    
    # Set axis limits
    ax.set_xlim(rmin - dr, rmax + dr)
    ax.set_ylim(-5, 185)
    
    # Add grid
    ax.grid(True, alpha=0.3, linestyle='--')
    
    # Add horizontal lines at common angles for reference
    ax.axhline(y=35.659, color='red', linestyle=':', alpha=0.7, label='θ = 35.659°')
    ax.axhline(y=60, color='orange', linestyle=':', alpha=0.7, label='θ = 60°')
    ax.axhline(y=109.47, color='red', linestyle=':', alpha=0.7, label='θ = 109.47°')
    ax.legend(loc='upper right')
    
    plt.tight_layout()
    
    if save_plot:
        plt.savefig(filename, dpi=300, bbox_inches='tight', 
                   facecolor='white', edgecolor='none')
        print(f"Radius-theta plot saved as {filename}")
    
    plt.show()

def plot_tbf_radius_theta_contour(tbf: Dict[str, Any], type_idx: int = 0, save_plot: bool = False, 
                                 filename: str = "tbf_radius_theta_contour.png") -> None:
    """
    Create a 2D contour plot showing radius-theta correlation where r1=r2=r.
    Theta is the angle between r1 and r2 vectors.
    
    Args:
        tbf: The TBF data structure returned by ReadTBF
        type_idx: Index of the type combination to plot (default: 0)
        save_plot: Whether to save the plot to a file (default: False)
        filename: Filename to save the plot (default: "tbf_radius_theta_contour.png")
    """
    # Extract parameters
    rmin = tbf['rmin']
    rmax = tbf['rmax']
    dr = tbf['dr']
    nr = tbf['nr']
    f3s = tbf['f3s'][type_idx]
    tmins = tbf['tmins']
    tmaxs = tbf['tmaxs']
    nt = tbf['nt']
    dt = tbf['dt']
    
    # Prepare data for radius-theta plotting
    r_values = []
    theta_values = []
    f3_values = []
    
    # Iterate through diagonal bins where r1 = r2 = r
    for ir in range(nr):
        # For r1 = r2, we use the diagonal where ir1 = ir2 = ir
        # MapIndex(ir, ir) gives us the bin index for r1 = r2 = r
        idx = MapIndex(ir, ir)
        r = rmin + dr * ir
        
        # Calculate theta bounds for this radius
        if nt[idx] > 0:
            for it in range(nt[idx]):
                theta = tmins[idx] + it * dt[idx]
                # Convert cosine to angle in degrees for better visualization
                angle_rad = math.acos(max(-1.0, min(1.0, theta)))
                angle_deg = math.degrees(angle_rad)
                
                # Only include valid angles and non-zero TBF values
                if 0 <= angle_deg <= 180 and f3s[idx][it] != 0:
                    r_values.append(r)
                    theta_values.append(angle_deg)
                    f3_values.append(f3s[idx][it])
    
    # Create contour plot
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Create a grid for contour plotting
    if len(r_values) > 0 and len(theta_values) > 0:
        # Create grid
        r_grid = np.linspace(min(r_values), max(r_values), 100)
        theta_grid = np.linspace(min(theta_values), max(theta_values), 100)
        R, Theta = np.meshgrid(r_grid, theta_grid)
        
        # Interpolate data onto the grid
        from scipy.interpolate import griddata
        points = np.column_stack((r_values, theta_values))
        
        # Create filled contour plot
        contourf = ax.contourf(R, Theta, griddata(points, f3_values, (R, Theta), method='cubic', fill_value=0), 
                              levels=50, cmap='YlGn', alpha=0.8)
        
        # Add contour lines
        contour_lines = ax.contour(R, Theta, griddata(points, f3_values, (R, Theta), method='cubic', fill_value=0), 
                                  levels=15, colors='black', alpha=0.3, linewidths=0.5)
        
        # Add contour labels
        ax.clabel(contour_lines, inline=True, fontsize=8, fmt='%.3f')
    
    else:
        # If no data, create empty plot with text
        ax.text(0.5, 0.5, 'No data available for contour plot', 
               transform=ax.transAxes, ha='center', va='center', fontsize=14)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
    
    # Customize the plot
    ax.set_xlabel('Radius r (Å)', fontsize=14, labelpad=10)
    ax.set_ylabel('Angle θ between r₁ and r₂ (degrees)', fontsize=14, labelpad=10)
    ax.set_title(f'Radius-Angle Contour Plot (r₁ = r₂ = r)', fontsize=16, pad=20)
    
    # Add colorbar
    cbar = plt.colorbar(contourf, ax=ax)
    cbar.set_label('TBF Value', fontsize=12)
    
    # Set axis limits
    if len(r_values) > 0:
        ax.set_xlim(min(r_values) - dr, max(r_values) + dr)
    else:
        ax.set_xlim(rmin - dr, rmax + dr)
    ax.set_ylim(-5, 185)
    
    # Add grid
    ax.grid(True, alpha=0.3, linestyle='--')
    
    # Add horizontal lines at common angles for reference
    ax.axhline(y=35.659, color='red', linestyle=':', alpha=0.7, label='θ = 35.659°')
    ax.axhline(y=60, color='orange', linestyle=':', alpha=0.7, label='θ = 60°')
    ax.axhline(y=109.47, color='red', linestyle=':', alpha=0.7, label='θ = 109.47°')
    ax.legend(loc='upper right')
    
    plt.tight_layout()
    
    if save_plot:
        plt.savefig(filename, dpi=300, bbox_inches='tight', 
                   facecolor='white', edgecolor='none')
        print(f"Radius-theta contour plot saved as {filename}")
    
    plt.show()

def plot_tbf_heatmap(tbf: Dict[str, Any], type_idx: int = 0, save_plot: bool = False, 
                    filename: str = "tbf_heatmap.png") -> None:
    """
    Create a 2D heatmap of the TBF data showing r1 vs r2 with color representing average TBF value.
    
    Args:
        tbf: The TBF data structure returned by ReadTBF
        type_idx: Index of the type combination to plot (default: 0)
        save_plot: Whether to save the plot to a file (default: False)
        filename: Filename to save the plot (default: "tbf_heatmap.png")
    """
    # Extract parameters
    rmin = tbf['rmin']
    rmax = tbf['rmax']
    dr = tbf['dr']
    nr = tbf['nr']
    f3s = tbf['f3s'][type_idx]
    nt = tbf['nt']
    
    # Create 2D arrays for plotting
    r1_grid = np.linspace(rmin, rmax, nr)
    r2_grid = np.linspace(rmin, rmax, nr)
    tbf_matrix = np.zeros((nr, nr))
    
    # Fill the matrix with average TBF values
    for ir1 in range(nr):
        for ir2 in range(nr):
            if ir1 >= ir2:  # Only consider r1 >= r2 due to symmetry
                idx = MapIndex(ir1, ir2)
                if nt[idx] > 0 and f3s[idx] is not None:
                    # Calculate average value for this bin
                    avg_value = sum(f3s[idx]) / len(f3s[idx])
                    tbf_matrix[ir1, ir2] = avg_value
                else:
                    tbf_matrix[ir1, ir2] = 0
            else:
                # For r1 < r2, use symmetry
                idx = MapIndex(ir2, ir1)
                if nt[idx] > 0 and f3s[idx] is not None:
                    avg_value = sum(f3s[idx]) / len(f3s[idx])
                    tbf_matrix[ir1, ir2] = avg_value
                else:
                    tbf_matrix[ir1, ir2] = 0
    
    # Create heatmap
    fig, ax = plt.subplots(figsize=(10, 8))
    im = ax.imshow(tbf_matrix, cmap='YlGn', origin='lower', 
                   extent=[rmin, rmax, rmin, rmax], aspect='auto')
    
    # Customize the plot
    ax.set_xlabel('r2 (Å)', fontsize=12)
    ax.set_ylabel('r1 (Å)', fontsize=12)
    ax.set_title(f'2D Heatmap\nAverage values over r3', fontsize=14)
    
    # Add colorbar
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('Average TBF Value', fontsize=12)
    
    # Add grid
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    if save_plot:
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        print(f"Heatmap saved as {filename}")
    
    plt.show()

# Example usage
if __name__ == "__main__":
    # Example of how to use the functions
    try:
        # Load TBF data
        tbf_data = ReadTBF("3bdf.dat")
        print_tbf_info(tbf_data)
        print("TBF data loaded successfully!")
        
        # Create 3D visualization with different plot types
        # print("\nCreating 3D scatter plot with auto threshold...")
        # plot_tbf_3d(tbf_data, type_idx=0, save_plot=True, filename="tbf_3d_scatter_auto.png", 
        #            plot_type='scatter', max_points=50000)
        
        # print("Creating 3D scatter plot with custom threshold...")
        # plot_tbf_3d(tbf_data, type_idx=0, save_plot=True, filename="tbf_3d_scatter_custom.png", 
        #            plot_type='scatter', max_points=50000, threshold=0.001)
        
        # print("Creating 3D surface plot...")
        # plot_tbf_3d(tbf_data, type_idx=0, save_plot=True, filename="tbf_3d_surface.png", 
        #            plot_type='surface')
        
        # print("Creating performance-optimized plot...")
        # plot_tbf_3d_performance(tbf_data, type_idx=0, save_plot=True, 
        #                       filename="tbf_3d_performance.png")
        
        # Create 2D radius-theta correlation plot
        print("Creating radius-theta correlation plot...")
        plot_tbf_radius_theta(tbf_data, type_idx=0, save_plot=True, filename="tbf_radius_theta.png")
        
        # Create 2D radius-theta contour plot
        print("Creating radius-theta contour plot...")
        plot_tbf_radius_theta_contour(tbf_data, type_idx=0, save_plot=True, filename="tbf_radius_theta_contour.png")
        
        # Create 2D heatmap
        print("Creating 2D heatmap...")
        plot_tbf_heatmap(tbf_data, type_idx=0, save_plot=True, filename="tbf_heatmap.png")
        
        print("All plots generated successfully!")
        # print("\nPlot types available:")
        # print("- scatter: Standard scatter plot (with point sampling for large datasets)")
        # print("- surface: Surface plot using interpolation (better for dense data)")
        # print("- contour: Contour plot (most efficient for very large datasets)")
        # print("- performance: Binned visualization for maximum performance")
        # print("- radius-theta: 2D plot showing correlation between radius and angle (r₁=r₂=r)")
        # print("\nThreshold options:")
        # print("- threshold=None: Auto-calculate as 1% of maximum absolute value")
        # print("- threshold=value: Use custom threshold value")
        # print("- Points with |TBF value| < threshold are filtered out")
        
    except FileNotFoundError:
        print("Input file '3bdf.dat' not found. Please provide the correct file path.")
    except Exception as e:
        print(f"Error: {e}")
        print("Make sure matplotlib and scipy are installed:")
        print("pip install matplotlib scipy")
