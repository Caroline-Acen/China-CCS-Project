"""
Figure 15 Supplementary: Historical Seismicity Map
China map showing historical earthquake events colored by magnitude.
"""

import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt

from config import DATA_DIR, NE_ADMIN0, NE_ADMIN1, chart_path, manuscript_font_bundle

HISTORICAL_EVENTS = DATA_DIR / "historical_events.xlsx"


def plot_historical_seismicity(
    historical_events,
    admin0_path,
    admin1_path,
    output_path="historical_events.png"
):
    """
    Create map showing historical seismicity events.

    Parameters:
        historical_events: Path to Excel file with latitude, longitude, magnitude columns
        admin0_path: Path to Natural Earth admin0 shapefile
        admin1_path: Path to Natural Earth admin1 shapefile
        output_path: Output image path
    """
    
    print("Loading historical seismicity data...")
    seismicity_df = pd.read_excel(historical_events)
    print(f"Loaded {len(seismicity_df)} events")
    print(f"Magnitude range: {seismicity_df['magnitude'].min():.1f} - {seismicity_df['magnitude'].max():.1f}")
    
    # Convert to GeoDataFrame
    seismicity_gdf = gpd.GeoDataFrame(
        seismicity_df,
        geometry=gpd.points_from_xy(seismicity_df["longitude"], seismicity_df["latitude"]),
        crs="EPSG:4326"
    )

    # Load boundaries
    print("Loading boundary data...")
    admin0 = gpd.read_file(admin0_path)
    admin1 = gpd.read_file(admin1_path)
    
    # Filter China
    china = admin0[admin0["ADMIN"] == "China"]
    
    if "admin" in admin1.columns:
        provinces = admin1[admin1["admin"] == "China"]
    elif "ADM0NAME" in admin1.columns:
        provinces = admin1[admin1["ADM0NAME"] == "China"]
    else:
        provinces = admin1[admin1["NAME"] == "China"]

    # Reproject
    target_crs = "ESRI:102012"
    china = china.to_crs(target_crs)
    provinces = provinces.to_crs(target_crs)
    seismicity_gdf = seismicity_gdf.to_crs(target_crs)

    # Create plot
    print("Creating plot...")
    fig, ax = plt.subplots(1, 1, figsize=(16, 12), dpi=300)
    mf = manuscript_font_bundle(fig.get_figwidth())
    
    # Plot boundaries
    china.boundary.plot(ax=ax, color="black", linewidth=1.5)
    provinces.boundary.plot(ax=ax, color="gray", linewidth=0.5, alpha=0.7)
    
    # Plot seismicity
    if not seismicity_gdf.empty:
        seismicity_size = seismicity_gdf['magnitude'] * 10
        
        scatter = ax.scatter(
            seismicity_gdf.geometry.x,
            seismicity_gdf.geometry.y,
            c=seismicity_gdf['magnitude'],
            s=seismicity_size,
            alpha=0.7,
            cmap='hot_r',
            edgecolors='black',
            linewidths=0.5
        )
        
        # Add colorbar
        cbar = plt.colorbar(scatter, ax=ax, shrink=0.8, pad=0.05)
        cbar.set_label('Earthquake Magnitude', fontsize=mf["colorbar_label"])
        cbar.ax.tick_params(labelsize=mf["colorbar_tick"])
    
    ax.axis('off')
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    
    print(f"[OK] Historical seismicity map saved to {output_path}")
    return fig


# Run function
if HISTORICAL_EVENTS.exists():
    plot_historical_seismicity(
        historical_events=str(HISTORICAL_EVENTS),
        admin0_path=str(NE_ADMIN0),
        admin1_path=str(NE_ADMIN1),
        output_path=chart_path("historical_seismicity_map.png"),
    )
else:
    print(f"[SKIP] Missing data file: {HISTORICAL_EVENTS}")
