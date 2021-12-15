from itertools import product

import cartopy.crs as ccrs
import numpy as np
import pandas as pd
import pylab as plt
import yaml
from cmcrameri import cm


def get_plotting_info_for_project(settings):
    def get_synth_stations(settings, wiggle=0.5):

        mode = settings["synth_stations_mode"]
        n = settings["synth_stations_n"]

        if mode == "grid":
            lons = np.linspace(-180, 180 - (360 / int(np.sqrt(n))), int(np.sqrt(n)))
            lats = np.linspace(-50, 50, int(np.sqrt(n)))
            station_locations = list(product(lons, lats))

        elif mode == "uniform":
            lons = np.random.uniform(low=0, high=180, size=n)
            lats = np.random.uniform(low=-50, high=50, size=n)
            station_locations = list(zip(lons, lats))

        elif mode == "partial_circle":
            n_total = settings["synth_stations_circle_max"]
            radius = settings["synth_stations_circle_radius"]
            n_used = settings["synth_stations_circle_n"]

            azimuths = np.linspace(0, 2 * np.pi, n_total)
            azimuths_used = azimuths[:n_used]

            lons = radius * np.cos(azimuths_used)
            lats = radius * np.sin(azimuths_used)
            station_locations = list(zip(lons, lats))

        elif mode == "file":

            df = pd.read_csv(
                "/data/cen/u254/sven/global_mfp/code/" + settings["synth_stations_file"]
            )
            lons = df["x"].values
            lats = df["y"].values
            station_locations = list(zip(lons, lats))

        station_locations = np.array(station_locations)

        return station_locations

    def generate_global_grid(settings):
        """
        Generate the grid geometry
        Returns coordinates of n grid_points as [[x0,y0,z0], [x1,y1,z1], ..., [xn,yn,zn]|
        """

        if settings["geometry_type"] == "cartesian":
            grid_limits_lon = settings["grid_limits_x"]
            grid_limits_lat = settings["grid_limits_y"]
        else:
            grid_limits_lon = settings["grid_limits_lon"]
            grid_limits_lat = settings["grid_limits_lat"]

        grid_spacing_in_deg = settings["grid_spacing"]

        n_gridpoints_lon = int(
            (grid_limits_lon[1] - grid_limits_lon[0]) / grid_spacing_in_deg
        )
        n_gridpoints_lat = int(
            (grid_limits_lat[1] - grid_limits_lat[0]) / grid_spacing_in_deg
        )

        # grid geometry
        grid_lon_coords = np.linspace(
            grid_limits_lon[0], grid_limits_lon[1], n_gridpoints_lon
        )
        grid_lat_coords = np.linspace(
            grid_limits_lat[0], grid_limits_lat[1], n_gridpoints_lat
        )
        grid_points = np.asarray(list(product(grid_lon_coords, grid_lat_coords)))

        return grid_points, grid_lon_coords, grid_lat_coords

    station_locations = get_synth_stations(settings)
    grid_points, grid_lon_coords, grid_lat_coords = generate_global_grid(
        settings=settings
    )
    source_loc = np.array(settings["synth_sources"]).T
    return grid_lon_coords, grid_lat_coords, station_locations, source_loc


def plot_beampowers_on_map(
    lons,
    lats,
    beampowers,
    settings,
    station_locations=None,
    outfile=None,
    source_loc=None,
    title=None,
    fig=None,
    ax=None,
    vmin=None,
    vmax=None,
    plot_station_locations=False,
    ax_id=None,
    LHS=False,
    eiwt_vel=None,
):

    xx, yy = np.meshgrid(lons, lats)
    trans = ccrs.PlateCarree()

    if settings["geometry_type"] == "geographic":
        ax.coastlines(resolution="10m", linewidth=0.5, color="#AAAAAA")

    if plot_station_locations:
        if settings["geometry_type"] == "geographic":
            ax.scatter(
                station_locations[:, 0],
                station_locations[:, 1],
                c="k",
                s=4,
                marker="^",
                lw=0,
                transform=trans,
                zorder=5,
            )
        else:
            ax.scatter(
                station_locations[:, 0],
                station_locations[:, 1],
                c="k",
                s=25,
                marker="^",
                lw=0,
                zorder=5,
                label="Station",
            )

    if settings["geometry_type"] == "geographic":
        if vmin is None and vmax is None:
            bp_absmax = np.abs(np.nanmax(beampowers))
            vmin = -bp_absmax
            vmax = bp_absmax
        pcm = ax.pcolormesh(
            xx,
            yy,
            beampowers.T,
            edgecolors="face",
            vmin=vmin,
            vmax=vmax,
            transform=trans,
            # cmap=cm.batlowW_r,
            cmap=cm.vik_r,
        )
    else:
        # force symmetric colorscale
        if vmin is None and vmax is None:
            bp_absmax = np.abs(np.nanmax(beampowers))
            vmin = -bp_absmax
            vmin = 0
            vmax = bp_absmax

        # compute max
        max_indices = np.unravel_index(
            np.argmax(beampowers, axis=None), beampowers.shape
        )
        lon_max = lons[max_indices[0]]
        lat_max = lats[max_indices[1]]
        # vmin = -.005
        # vmax = .005
        pcm = ax.pcolormesh(
            xx,
            yy,
            beampowers.T,
            edgecolors="face",
            vmin=vmin,
            vmax=vmax,
            # cmap=cm.batlowW_r,
            cmap=cm.vik_r,
            shading="nearest",
        )
        # pcm = ax.contourf(xx, yy, beampowers.T, levels=np.linspace(-1, 1, 50), vmin=vmin, vmax=vmax, cmap=roma_map) # , norm=LogNorm())
        ax.set_xticks([-50, 0, 50])
        ax.set_xticklabels([-50, 0, 50], fontsize=6)
        ax.set_yticks([-50, 0, 50])
        ax.set_yticklabels([-50, 0, 50], fontsize=6)
        if ax_id == 0:
            if LHS:
                ax.set_ylabel("Distance [km]", fontsize=6)
                ax.set_xticklabels([])
                ax.set_yticklabels(["-50\n50", 0, 50], fontsize=6)
            else:
                ax.set_xticklabels([])
                ax.set_yticklabels(["-50\n50", 0, 50], fontsize=6)
        elif ax_id == 1:
            ax.set_xticklabels([])
            ax.set_yticklabels([])
        elif ax_id == 2:
            if LHS:
                ax.set_ylabel("Distance [km]", fontsize=6)
                ax.set_yticklabels([-50, 0, ""], fontsize=6)
                ax.set_xticklabels([-50, 0, "50/-50"], fontsize=6)
                # ax.set_xlabel("Distance [km]", fontsize=6)
            else:
                ax.set_yticklabels([-50, 0, ""], fontsize=6)
                ax.set_xticklabels([-50, 0, "50/-50"], fontsize=6)
                # ax.set_xlabel("Distance [km]", fontsize=6)
        elif ax_id == 3:
            ax.set_yticklabels([])
            # ax.set_xlabel("Distance [km]", fontsize=6)
            ax.set_xticklabels(["", 0, 50], fontsize=6)

        ax.set_xlim(-50, 50)
        ax.set_ylim(-50, 50)

        if ax_id is None:
            ax.set_xlabel("Distance [km]", fontsize=6)
            if LHS:
                ax.set_ylabel("Distance [km]", fontsize=6)
        dist_to_true_src = np.linalg.norm(
            [lon_max - source_loc[0], lat_max - source_loc[1]]
        )
        if eiwt_vel:
            t = ax.text(
                0.07,
                0.93,
                f"v={eiwt_vel}km/s\n"
                + r"$\Delta \vec{x}$ = "
                + f"{dist_to_true_src:.1f}km",
                ha="left",
                va="top",
                fontsize=6,
                # fontweight="bold",
                transform=ax.transAxes,
            )
            t.set_bbox(dict(facecolor="w", alpha=0.50, edgecolor="None"))
        else:
            t = ax.text(
                0.05,
                0.95,
                f"Our approach\n"
                + r"$\Delta \vec{x}$ = "
                + f"{dist_to_true_src:.1f}km",
                ha="left",
                va="top",
                fontsize=6,
                # fontweight="bold",
                transform=ax.transAxes,
            )
            t.set_bbox(dict(facecolor="w", alpha=0.50, edgecolor="None"))

        if LHS and ax_id == 0:
            ax.text(
                0, 1.05, "a)", ha="center", va="bottom", transform=ax.transAxes,
            )
            ax.text(
                1, 1.05, "broadband", ha="center", va="bottom", transform=ax.transAxes,
            )
        if ax_id == 0 and not LHS:
            ax.text(
                0, 1.05, "c)", ha="center", va="bottom", transform=ax.transAxes,
            )
            ax.text(
                1,
                1.05,
                "narrowband (0.13-0.15Hz)",
                ha="center",
                va="bottom",
                transform=ax.transAxes,
            )
        if ax_id is None and LHS:
            ax.text(
                0, 1.025, "b)", ha="center", va="bottom", transform=ax.transAxes,
            )
        if ax_id is None and not LHS:
            ax.text(
                0, 1.025, "d)", ha="center", va="bottom", transform=ax.transAxes,
            )
        # pcm = ax.pcolormesh(xx, yy, beampowers.T, edgecolors='face', vmin=-bp_absmax, vmax=bp_absmax, norm=SymLogNorm(linthresh=1E-3), cmap=roma_map)
    # ax.stock_img()

    # x0, y0, w, h = ax.get_position().bounds
    # cb_ax = fig.add_axes([x0 + w + 0.025 * w, y0, 0.025 * w, h])
    # plt.colorbar(pcm, cax=cb_ax)

    # plot max
    if settings["do_synth"]:

        # tdist = ax.text(
        #     source_loc[0] + 5,
        #     source_loc[1] + 5,
        #     r"$\Delta \vec{x}$ = " + f"{dist_to_true_src:.2f}km",
        #     fontsize=6,
        #     va="bottom",
        #     ha="left",
        # )
        if settings["geometry_type"] == "geographic":
            ax.scatter(
                lon_max,
                lat_max,
                facecolors="#E2001A",
                edgecolors="w",
                linewidth=0.5,
                s=20,
                marker="o",
                label="Beampower Peak",
                transform=trans,
            )
        elif settings["geometry_type"] == "cartesian":
            ax.scatter(
                lon_max,
                lat_max,
                facecolors="#E2001A",
                edgecolors="w",
                linewidth=0.5,
                s=20,
                marker="o",
                label="Beampower Peak",
                zorder=20,
            )

    if source_loc is not None:
        if settings["geometry_type"] == "geographic":
            ax.scatter(
                source_loc[0],
                source_loc[1],
                edgecolors="k",
                c="magenta",
                linewidth=0.5,
                s=50,
                marker="*",
                transform=trans,
            )
        else:
            ax.scatter(
                source_loc[0, :],
                source_loc[1, :],
                edgecolors="k",
                # c="#e2001a",
                c="w",
                linewidth=0.5,
                s=100,
                marker="*",
                label="Synth. Source",
            )

    if title and settings["geometry_type"] == "geographic":
        ax.set_title(title)

    # ax.legend()

    # from cartopy.io import shapereader
    # ax.add_geometries(list(shapereader.Reader('/home/zmaw/u254070/.local/share/cartopy/shapefiles/natural_earth/physical/ne_10m_coastline.shp').geometries()), crs_use)

    # import cartopy.feature as cfeat
    # ax.add_feature(cfeat.GSHHSFeature)

    # ax.gridlines(crs=crs_use)

    # custom lim to check results
    # chile lat -35.47, lon -72.91
    # ax.set_extent([-76, -72, -40, -25], crs=trans)
    # ax.set_ylim(-40, -30)

    if ax_id is None and LHS:
        lgnd = ax.legend(loc="upper right", fontsize=6)

    if outfile:
        fig.savefig(outfile, dpi=300, transparent=False)
    else:
        plt.show()
    plt.close(fig)

    return pcm


plt.rcParams.update({"font.size": 8})

with open(
    "/data/cen/u254/sven/global_mfp/projects/173_redevelop_synth_test_larger_array_GF/1635927790_settings.yml",
    "r",
) as stream:
    settings = yaml.safe_load(stream)

# fig, axs = plt.subplots(2, 2)
fig = plt.figure(figsize=(5, 5))
fig.subplots_adjust(left=0.09, right=0.915, top=0.95)
outer_grid = fig.add_gridspec(2, 2)
# left plots (broadband)
# e-iwt - subplots for vconst=2.5, 2.6, 2.7, 2.8
inner_grid_left = outer_grid[0, 0].subgridspec(2, 2, wspace=0, hspace=0)
axs_eiwt_broadband = inner_grid_left.subplots()
# eiwt_vels = [2.5, 2.6, 2.7, 2.8]
eiwt_vels = [3.0, 3.2, 3.4, 3.6]
for ax_id, (ax, eiwt_vel) in enumerate(zip(axs_eiwt_broadband.flatten(), eiwt_vels)):
    bp_fn = f"/data/cen/u254/sven/global_mfp/projects/173_redevelop_synth_test_larger_array_v{eiwt_vel}/out/out_1548982800.0_3880_None_0_0_False.npy"
    # bp_fn = f"/data/cen/u254/sven/global_mfp/projects/172_redo_synth_fig2_vconst{eiwt_vel}/out/out_1217356935.0_3880_None_0_0_False.npy"
    beampowers = np.load(bp_fn)
    beampowers /= np.max(np.abs(beampowers))
    print(ax)
    (
        grid_lon_coords,
        grid_lat_coords,
        station_locations,
        source_loc,
    ) = get_plotting_info_for_project(settings)
    pcm = plot_beampowers_on_map(
        lons=grid_lon_coords,
        lats=grid_lat_coords,
        beampowers=beampowers,
        settings=settings,
        station_locations=station_locations,
        source_loc=source_loc,
        ax=ax,
        vmin=-1,
        vmax=1,
        plot_station_locations=True,
        ax_id=ax_id,
        LHS=True,
        eiwt_vel=eiwt_vel,
    )
# GF
# bp_fn = "/data/cen/u254/sven/global_mfp/projects/159_synth_test_GF/out/out_1548982800.0_3880_None_0_0.npy"
with open(
    "/data/cen/u254/sven/global_mfp/projects/173_redevelop_synth_test_larger_array_GF_densergrid/1635941547_settings.yml",
    "r",
) as stream:
    settings = yaml.safe_load(stream)
bp_fn = f"/data/cen/u254/sven/global_mfp/projects/173_redevelop_synth_test_larger_array_GF_densergrid/out/out_1548982800.0_3880_[0.01, 0.5]_0_0_False.npy"
beampowers = np.load(bp_fn)
beampowers /= np.max(np.abs(beampowers))
ax_GF_broadband = fig.add_subplot(outer_grid[1, 0])
(
    grid_lon_coords,
    grid_lat_coords,
    station_locations,
    source_loc,
) = get_plotting_info_for_project(settings)
pcm = plot_beampowers_on_map(
    lons=grid_lon_coords,
    lats=grid_lat_coords,
    beampowers=beampowers,
    settings=settings,
    station_locations=station_locations,
    source_loc=source_loc,
    ax=ax_GF_broadband,
    vmin=-1,
    vmax=1,
    plot_station_locations=True,
    LHS=True,
)

# right plots (narrowband)
# e-iwt
with open(
    "/data/cen/u254/sven/global_mfp/projects/173_redevelop_synth_test_larger_array_GF/1635927790_settings.yml",
    "r",
) as stream:
    settings = yaml.safe_load(stream)
inner_grid_left = outer_grid[0, 1].subgridspec(2, 2, wspace=0, hspace=0)
axs_eiwt_narrowband = inner_grid_left.subplots()
for ax_id, (ax, eiwt_vel) in enumerate(zip(axs_eiwt_narrowband.flatten(), eiwt_vels)):
    # bp_fn = f"/data/cen/u254/sven/global_mfp/projects/159_synth_test_vconst{eiwt_vel}/out/out_1548982800.0_3880_[0.13, 0.15]_0_0.npy"
    # bp_fn = f"/data/cen/u254/sven/global_mfp/projects/172_redo_synth_fig2_vconst{eiwt_vel}/out/out_1217356935.0_3880_[0.13, 0.15]_0_0_False.npy"
    bp_fn = f"/data/cen/u254/sven/global_mfp/projects/173_redevelop_synth_test_larger_array_v{eiwt_vel}/out/out_1548982800.0_3880_[0.13, 0.15]_0_0_False.npy"
    beampowers = np.load(bp_fn)
    beampowers /= np.max(np.abs(beampowers))
    (
        grid_lon_coords,
        grid_lat_coords,
        station_locations,
        source_loc,
    ) = get_plotting_info_for_project(settings)
    pcm = plot_beampowers_on_map(
        lons=grid_lon_coords,
        lats=grid_lat_coords,
        beampowers=beampowers,
        settings=settings,
        station_locations=station_locations,
        source_loc=source_loc,
        ax=ax,
        vmin=-1,
        vmax=1,
        plot_station_locations=True,
        ax_id=ax_id,
        LHS=False,
        eiwt_vel=eiwt_vel,
    )
# GF
with open(
    "/data/cen/u254/sven/global_mfp/projects/173_redevelop_synth_test_larger_array_GF_densergrid/1635941547_settings.yml",
    "r",
) as stream:
    settings = yaml.safe_load(stream)
bp_fn = f"/data/cen/u254/sven/global_mfp/projects/173_redevelop_synth_test_larger_array_GF_densergrid/out/out_1548982800.0_3880_[0.13, 0.15]_0_0_False.npy"
beampowers = np.load(bp_fn)
beampowers /= np.max(np.abs(beampowers))
ax_GF_narrowband = fig.add_subplot(outer_grid[1, 1])
(
    grid_lon_coords,
    grid_lat_coords,
    station_locations,
    source_loc,
) = get_plotting_info_for_project(settings)
pcm = plot_beampowers_on_map(
    lons=grid_lon_coords,
    lats=grid_lat_coords,
    beampowers=beampowers,
    settings=settings,
    station_locations=station_locations,
    source_loc=source_loc,
    ax=ax_GF_narrowband,
    vmin=-1,
    vmax=1,
    plot_station_locations=True,
    LHS=False,
)

cb_ax = fig.add_axes([0.915 + 0.015, 0.3, 0.015, 0.4])
cbar = plt.colorbar(pcm, cax=cb_ax)  # , extend="min")
cbar.set_ticks([-1, 0, 1])
cbar.set_label(label="Normalised Beampower", fontsize=6, labelpad=-3)
cbar.ax.tick_params(labelsize=6)

fig.savefig("/data/cen/u254/sven/global_mfp/figure_scripts/png/fig2.png", dpi=300)
