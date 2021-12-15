# only plot stack, if something to stack
# from sven_utils import suColors
# my_cmap = suColors.get_custom_cmap(name='devon', inverted=True)
# if len(beampowers_per_start_time) > 1:

from cmcrameri import cm


def get_plotting_info_for_project(settings):
    def get_synth_stations(settings, wiggle=0.5):
        from itertools import product

        import numpy as np

        mode = settings["synth_stations_mode"]
        n = settings["synth_stations_n"]

        if mode == "grid":
            lons = np.linspace(-180, 180 - (360 / int(np.sqrt(n))), int(np.sqrt(n)))
            lats = np.linspace(-75, 50, int(np.sqrt(n)))
            station_locations = list(product(lons, lats))

        elif mode == "uniform":
            lons = np.random.uniform(low=0, high=180, size=n)
            lats = np.random.uniform(low=-75, high=75, size=n)
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
            import pandas as pd

            df = pd.read_csv(
                "/data/cen/u254/sven/global_mfp/code/" + settings["synth_stations_file"]
            )
            lons = df["x"].values
            lats = df["y"].values
            station_locations = list(zip(lons, lats))

        # if wiggle != 0:
        # station_locations = [[sta_lon+np.random.uniform(-wiggle, wiggle), sta_lat+np.random.uniform(-wiggle, wiggle)] for sta_lon, sta_lat in product(lons, lats)]

        station_locations = np.array(station_locations)

        return station_locations

    def generate_global_grid(settings):
        """
        Generate the grid geometry
        Returns coordinates of n grid_points as [[x0,y0,z0], [x1,y1,z1], ..., [xn,yn,zn]|
        """
        from itertools import product

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

        # exclude grid points on land?
        # from global_land_mask import globe
        # for gp in grid_points:
        #     if globe.is_land(gp[1], gp[0])

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
    nonorm=False,
):
    import cartopy.crs as ccrs
    import numpy as np
    import pylab as plt

    # plt.style.use("dracula")
    # crs_use = ccrs.Robinson()
    # convert lons, lats to corresponding CRS

    xx, yy = np.meshgrid(lons, lats)

    # ax = plt.axes(projection=ccrs.Orthographic(central_longitude=source_loc[0], central_latitude=source_loc[1]))
    # _map = Basemap(projection='eck4',lon_0=0,resolution='c', ax=ax)
    # _map.drawcoastlines(linewidth=.5, color='k')
    # _map.drawparallels(np.arange(-90.,120.,30.))
    # _map.drawmeridians(np.arange(0.,420.,60.))
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
                ec="w",
                s=20,
                marker="^",
                lw=0.5,
                zorder=5,
                label="Station",
            )

    # xx, yy = _map(xx_mg, yy_mg)
    from matplotlib.colors import LogNorm, SymLogNorm

    if settings["geometry_type"] == "geographic":
        # bp_absmax = np.abs(np.nanmax(beampowers))
        # bp_absmin = np.abs(np.nanmin(beampowers))
        if vmin is None and vmax is None:
            bp_absmax = np.abs(np.nanmax(beampowers))
            vmin = -bp_absmax
            vmax = bp_absmax
        # pcm = ax.pcolormesh(xx, yy, beampowers.T, edgecolors='face', vmin=-bp_absmax, vmax=bp_absmax, transform=trans, norm=SymLogNorm(linthresh=1E-3), cmap='RdBu')
        # pcm = ax.pcolormesh(xx, yy, beampowers.T, edgecolors='face', vmin=bp_absmin, vmax=bp_absmax, transform=trans, norm=LogNorm(), cmap='magma')
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
        # vmin = -.005
        # vmax = .005
        if nonorm:
            pcm = ax.pcolormesh(
                xx,
                yy,
                beampowers.T,
                edgecolors="face",
                vmin=1e-4,
                vmax=vmax,
                # cmap=cm.batlowW_r,
                cmap=cm.vik_r,
                shading="nearest",
                norm=LogNorm(),
            )
            x0, y0, w, h = ax.get_position().bounds
            cb_ax = fig.add_axes([x0 + 0.9 * w, y0 + 0.2 * h, w * 0.05, h * 0.6])
            cbar = plt.colorbar(pcm, cax=cb_ax)
            cbar.set_ticks([-1, 0, 1])
            cbar.set_label(label="Normalised Beampower", fontsize=6)
            cbar.ax.tick_params(labelsize=6)
        else:
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
            )  # , norm=LogNorm())
        # pcm = ax.contourf(xx, yy, beampowers.T, levels=np.linspace(-1, 1, 50), vmin=vmin, vmax=vmax, cmap=roma_map) # , norm=LogNorm())
        ax.set_xticks([-50, 0, 50])
        ax.set_xticklabels([-50, 0, 50], fontsize=6)
        ax.set_yticks([-50, 0, 50])
        ax.set_yticklabels([-50, 0, 50], fontsize=6)
        figlbl = ""
        if ax_id == 0:
            figlbl = "a)"
            labeltext = "Station density\n(x4 on right side)"
            ax.set_xlabel("Distance [km]", fontsize=6)
            ax.set_ylabel("Distance [km]", fontsize=6)
            ax.scatter(
                station_locations[3:, 0],
                station_locations[3:, 1],
                c="k",
                ec="w",
                s=60,
                marker="^",
                lw=0.75,
                zorder=5,
            )
            # if LHS:
            #     ax.set_ylabel("Distance [km]", fontsize=6)
            #     ax.set_xticklabels([])
            #     ax.set_yticklabels([-50, 0, 50], fontsize=6)
            # else:
            #     ax.set_xticklabels([])
            #     ax.set_yticklabels([-50, 0, 50], fontsize=6)
        elif ax_id == 1:
            figlbl = "a)"
            labeltext = "Right: quadruple station density"
            ax.set_xticklabels([])
            ax.set_yticklabels([])
        elif ax_id == 2:
            figlbl = "b)"
            labeltext = "Two sources (narrowband)"
            ax.set_yticklabels([])
            # ax.set_yticks([])
            ax.set_xlabel("Distance [km]", fontsize=6)
            ax.set_xticklabels([-50, 0, 50], fontsize=6)
            # ax.set_xlabel("Distance [km]", fontsize=6)

        elif ax_id == 3:
            figlbl = "c)"
            labeltext = "Two sources (broadband)"
            ax.set_yticklabels([])
            # ax.set_yticks([])
            ax.set_xlabel("Distance [km]", fontsize=6)
            ax.set_xticklabels([-50, 0, 50], fontsize=6)
            from scipy.ndimage import maximum_filter

            maxima = beampowers == maximum_filter(beampowers, size=50, mode="constant")
            localmaxlons = lons[np.where(maxima)[0]]
            localmaxlats = lats[np.where(maxima)[1]]
            maxima_values = beampowers[maxima]
            localmaxlons[maxima_values <= 0.5] = np.nan
            localmaxlats[maxima_values <= 0.5] = np.nan
            # localmaxlons[localmaxlats <= -80] = np.nan
            # localmaxlats[localmaxlats <= -80] = np.nan
            # localmaxlats[localmaxlons <= -5] = np.nan
            # localmaxlons[localmaxlons <= -5] = np.nan
            # localmaxlats[localmaxlons >= 5] = np.nan
            # localmaxlons[localmaxlons >= 5] = np.nan
            ax.scatter(
                localmaxlons,
                localmaxlats,
                facecolors="#E2001A",
                edgecolors="#E2001A",
                linewidth=0.5,
                s=5,
                marker="o",
                label="Beampower Peak",
                alpha=1,
                zorder=10,
            )
        ax.set_xlim(-50, 50)
        ax.set_ylim(-50, 50)
        ltextx, ltexty = 0.05, 0.95
        if ax_id is None:
            ax.set_xlabel("Distance [km]", fontsize=6)
            if LHS:
                ax.set_ylabel("Distance [km]", fontsize=6)
        if ax_id == 100:
            labeltext = "zoom-in"
            ltextx += 0.02
            ltexty -= 0.02
        t = ax.text(
            ltextx,
            ltexty,
            labeltext,
            ha="left",
            va="top",
            fontsize=6,
            # fontweight="bold",
            transform=ax.transAxes,
        )
        t.set_bbox(dict(facecolor="w", alpha=0.75, edgecolor="None"))
        ax.text(
            0.00,
            1.01,
            figlbl,
            ha="left",
            va="bottom",
            fontsize=10,
            # fontweight="bold",
            transform=ax.transAxes,
            zorder=100,
        )
        # if LHS and ax_id == 0:
        #     ax.text(
        #         0, 1.05, "a)", ha="center", va="bottom", transform=ax.transAxes,
        #     )
        #     # ax.text(
        #     #     1, 1.05, "broadband", ha="center", va="bottom", transform=ax.transAxes,
        #     # )

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
        max_indices = np.unravel_index(
            np.argmax(beampowers, axis=None), beampowers.shape
        )
        lon_max = lons[max_indices[0]]
        lat_max = lats[max_indices[1]]
        if settings["geometry_type"] == "geographic":
            ax.scatter(
                lon_max,
                lat_max,
                facecolors="#E2001A",
                edgecolors="w",
                linewidth=0.5,
                s=50,
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
                zorder=10,
            )

    if source_loc is not None:
        if settings["geometry_type"] == "geographic":
            ax.scatter(
                source_loc[0],
                source_loc[1],
                edgecolors="k",
                c="magenta",
                linewidth=0.75,
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

    if ax_id == 100:
        # 65,-66.33
        ax.set_xlim(65 - 2.5, 65 + 2.5)
        ax.set_ylim(-66.33 - 2.5, -66.33 + 2.5)
        ax.set_xticks([64, 66])
        ax.set_yticks([-64, -66, -68])
        ax.set_xticklabels([64, 66], fontsize=6)
        ax.set_yticklabels([-64, -66, -68], fontsize=6)
        # 65 - 2.5, 65 + 2.5
        # -66.33 - 2.5, -66.33 + 2.5
        # 65,-66.33

    if ax_id == 0:
        # lgnd = ax.legend(loc="upper right", fontsize=6)
        pass

    if outfile:
        fig.savefig(outfile, dpi=300, transparent=False)
    else:
        plt.show()
    plt.close(fig)

    return pcm


import numpy as np
import pylab as plt
import yaml

plt.rcParams.update({"font.size": 8})

with open(
    "/data/cen/u254/sven/global_mfp/projects/173_redo_quad/1636027988_settings.yml",
    "r",
) as stream:
    settings = yaml.safe_load(stream)

# fig, axs = plt.subplots(2, 2)
fig = plt.figure(figsize=(6, 2.25))
fig.subplots_adjust(left=0.075, bottom=0.075, right=0.915, top=0.95)
outer_grid = fig.add_gridspec(1, 3, wspace=0.1, hspace=0.1)
# left plots (broadband)
# e-iwt - subplots for vconst=2.5, 2.6, 2.7, 2.8
# inner_grid_left = outer_grid[0, 0].subgridspec(2, 2, wspace=0, hspace=0)
# axs_eiwt_broadband = inner_grid_left.subplots()
# ax_notreatment = fig.add_subplot(outer_grid[0])
# # bp_fn = f"/data/cen/u254/sven/global_mfp/projects/160_ampl_treatment_none/out/out_1548982800.0_3880_[0.13, 0.15]_0_0.npy"
# # bp_fn = f"/data/cen/u254/sven/global_mfp/projects/167_synth/out/out_1548982800.0_3880_[0.13, 0.15]_0_0_False.npy"
# bp_fn = f"/data/cen/u254/sven/global_mfp/projects/170_redo_synths_fig1_broadband_GF/out/out_1548982800.0_3880_[0.13, 0.15]_0_0_False.npy"
# # bp_fn = f"/data/cen/u254/sven/global_mfp/projects/160_ampl_treatment_spreading_attenuation_7_corciulo/out/out_1548982800.0_3880_[0.13, 0.15]_0_0.npy"
# beampowers = np.load(bp_fn)
# beampowers /= np.max(np.abs(beampowers))
# (
#     grid_lon_coords,
#     grid_lat_coords,
#     station_locations,
#     source_loc,
# ) = get_plotting_info_for_project(settings)
# pcm = plot_beampowers_on_map(
#     lons=grid_lon_coords,
#     lats=grid_lat_coords,
#     beampowers=beampowers,
#     settings=settings,
#     station_locations=station_locations,
#     source_loc=source_loc,
#     ax=ax_notreatment,
#     vmin=-1,
#     vmax=1,
#     plot_station_locations=True,
#     ax_id=0,
#     LHS=True,
#     # nonorm=True,
#     # fig=fig,
# )

# ax_notreatment.plot(
#     [0.335, 0.335 + 0.09], [0.675, 0.61], transform=fig.transFigure, c="k", lw=0.5,
# )

# 65 - 2.5, 65 + 2.5
# -66.33 - 2.5, -66.33 + 2.5
# ax_notreatment.plot(
#     [65 - 2.5, 65 - 2.5, 65 + 2.5, 65 + 2.5, 65 - 2.5],
#     [-66.33 - 2.5, -66.33 + 2.5, -66.33 + 2.5, -66.33 - 2.5, -66.33 - 2.5],
#     c="#E2001A",
#     lw=1,
#     zorder=100,
# )

# ax_notreatment_zoom.set_xlim()

# time-domain top-right
ax_timedomain = fig.add_subplot(outer_grid[0])
ax_timedomain.set_aspect("equal")
bp_fn = f"/data/cen/u254/sven/global_mfp/projects/173_redo_quad/out/out_1548982800.0_3880_[0.13, 0.15]_0_0_False.npy"
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
    ax=ax_timedomain,
    vmin=-1,
    vmax=1,
    plot_station_locations=True,
    ax_id=0,
    LHS=False,
)


with open(
    "/data/cen/u254/sven/global_mfp/projects/173_redo_two_source_12/1636104599_settings.yml",
    "r",
) as stream:
    settings = yaml.safe_load(stream)

(
    grid_lon_coords,
    grid_lat_coords,
    station_locations,
    source_loc,
) = get_plotting_info_for_project(settings)

# phase-corr bottom-left
ax_phasecorr = fig.add_subplot(outer_grid[1])
ax_phasecorr.set_aspect("equal")
bp_fn = f"/data/cen/u254/sven/global_mfp/projects/173_redo_two_source_12/out/out_1548982800.0_3880_[0.13, 0.15]_0_0_False.npy"
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
    ax=ax_phasecorr,
    vmin=-1,
    vmax=1,
    plot_station_locations=True,
    ax_id=2,
    LHS=True,
)

# whitening bottom-right
ax_whitening = fig.add_subplot(outer_grid[2])
ax_whitening.set_aspect("equal")
bp_fn = f"/data/cen/u254/sven/global_mfp/projects/173_redo_two_source_12/out/out_1548982800.0_3880_[0.01, 0.5]_0_0_False.npy"
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
    ax=ax_whitening,
    vmin=-1,
    vmax=1,
    plot_station_locations=True,
    ax_id=3,
    LHS=False,
)


x0, y0, w, h = ax_whitening.get_position().bounds
#  cb_ax = fig.add_axes([0.915 + 0.015, 0.3, 0.015, 0.4])
cb_ax = fig.add_axes([x0 + w + 0.05 * w, y0, 0.05 * w, h])
cbar = plt.colorbar(pcm, cax=cb_ax) # , extend="min")
cbar.set_ticks([-1, 0, 1])
cbar.set_label(label="Normalised Beampower", fontsize=6)
cbar.ax.tick_params(labelsize=6)

fig.savefig("/data/cen/u254/sven/global_mfp/figure_scripts/png/fig4.png", dpi=300)
