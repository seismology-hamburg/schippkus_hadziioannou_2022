# plot earthquake real data example

# only plot stack, if something to stack
# from sven_utils import suColors
# my_cmap = suColors.get_custom_cmap(name='devon', inverted=True)
# if len(beampowers_per_start_time) > 1:

import cartopy.crs as ccrs
import cartopy.feature
import obspy
from cmcrameri import cm
from matplotlib.pyplot import axis


def get_plotting_info_for_project(settings):
    import obspy

    def get_station_locs_from_staxml(st, settings):
        import obspy

        station_locations = []
        for tr in st:
            inv = obspy.read_inventory(
                f"./data/stations/{tr.stats.network}.{tr.stats.station}.xml",
                format="STATIONXML",
            )
            lat, lon = inv[0][0][0].latitude, inv[0][0][0].longitude
            station_locations.append([lon, lat])
        return np.array(station_locations)

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

    #
    st = obspy.read("./data/deconv/CI/LHZ_chino.mseed")
    st.merge(fill_value=0)
    station_locations = get_station_locs_from_staxml(st, settings)
    # print(station_locations)
    grid_points, grid_lon_coords, grid_lat_coords = generate_global_grid(
        settings=settings
    )
    # source_loc = np.array(settings["synth_sources"]).T
    source_loc = [-117.766, 33.949]
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
    import cartopy.crs as ccrs
    import numpy as np
    import pylab as plt

    # plt.style.use("dracula")
    # crs_use = ccrs.Robinson()
    # convert lons, lats to corresponding CRS

    xx, yy = np.meshgrid(lons, lats)

    max_indices = np.unravel_index(np.argmax(beampowers, axis=None), beampowers.shape)
    lon_max = lons[max_indices[0]]
    lat_max = lats[max_indices[1]]
    dist = (
        obspy.geodetics.gps2dist_azimuth(
            lat1=source_loc[1], lon1=source_loc[0], lat2=lat_max, lon2=lon_max
        )[0]
        / 1000
    )

    # ax = plt.axes(projection=ccrs.Orthographic(central_longitude=source_loc[0], central_latitude=source_loc[1]))
    # _map = Basemap(projection='eck4',lon_0=0,resolution='c', ax=ax)
    # _map.drawcoastlines(linewidth=.5, color='k')
    # _map.drawparallels(np.arange(-90.,120.,30.))
    # _map.drawmeridians(np.arange(0.,420.,60.))
    trans = ccrs.PlateCarree()

    if settings["geometry_type"] == "geographic":
        ax.coastlines(resolution="10m", linewidth=0.5, color="k")

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
        # cnt = ax.contour(
        #     xx,
        #     yy,
        #     beampowers.T,
        #     levels=[0.95],
        #     transform=trans,
        #     colors="#E2001A",
        #     linewidths=1,
        # )
        # ax.set_xlim(
        #     settings["grid_limits_lon"][0] + 0.75, settings["grid_limits_lon"][1] - 1
        # )
        # ax.set_ylim(settings["grid_limits_lat"][0], settings["grid_limits_lat"][1])
        # plot_lims = (
        #     np.round(source_loc[0] - 2, 0),
        #     np.round(source_loc[0] + 2, 0),
        #     np.round(source_loc[1] - 1.5, 0),
        #     np.round(source_loc[1] + 1.5, 0),
        # )
        plot_lims = (-119, -116.5, 33, 35)
        ax.set_xlim(*plot_lims[:2])
        ax.set_ylim(*plot_lims[2:])

        if ax_id == 0:
            ax.set_xticks([-119, -118, -117])
            ax.set_yticks([33, 34, 35])
            ax.set_xticklabels(
                [
                    r"119$^\circ$W",
                    r"118$^\circ$W",
                    r"117$^\circ$W",
                    # r"116$^\circ$W/119$^\circ$W",
                ],
                fontsize=6,
            )
            ax.set_yticklabels(
                [r"33$^\circ$N", r"34$^\circ$N", r"35$^\circ$N"], fontsize=6
            )

            # overview_ax
            x0, y0, w, h = ax.get_position().bounds
            overview_ax = fig.add_axes(
                [x0 + 0.65 * w, y0 + 0.63 * h, 0.35 * w, 0.35 * h],
                projection=ccrs.NearsidePerspective(
                    central_latitude=33.949,
                    central_longitude=-117.766,
                    satellite_height=500_000,
                ),
            )
            # overview_ax.set_frame_on(False)

            overview_ax.plot(
                [plot_lims[0], plot_lims[0], plot_lims[1], plot_lims[1], plot_lims[0]],
                [plot_lims[2], plot_lims[3], plot_lims[3], plot_lims[2], plot_lims[2]],
                c="k",
                lw=0.75,
                transform=trans,
            )

            overview_ax.add_feature(cartopy.feature.LAND, color="w")
            overview_ax.add_feature(cartopy.feature.OCEAN, color="w")
            overview_ax.add_feature(cartopy.feature.STATES, lw=0.1, edgecolor="#CCCCCC")
            overview_ax.add_feature(cartopy.feature.BORDERS, lw=0.25, color="#AAAAAA")
            overview_ax.set_global()

            overview_ax.coastlines(lw=0.25)
            # overview_ax.scatter(
            #     33.949,
            #     -117.766,
            #     edgecolors="k",
            #     c="w",
            #     linewidth=0.5,
            #     s=25,
            #     marker="*",
            #     transform=trans,
            #     zorder=5,
            # )

            label_text = (
                "Standard approach\n" + r"$\Delta \vec{x}$ = " + f"{dist:0.1f}km"
            )

        # ax.add_feature(cartopy.feature.STATES, lw=0.1, edgecolor="#CCCCCC")
        ax.add_feature(cartopy.feature.BORDERS, lw=0.5, color="k")

        if ax_id == 1:
            label_text = (
                "Our approach - Explosion source\n"
                + r"$\Delta \vec{x} = $"
                + f"{dist:0.1f}km"
            )
            ax.set_xticks([-119, -118, -117])
            ax.set_xticklabels(
                [r"119$^\circ$W", r"118$^\circ$W", r"117$^\circ$W"], fontsize=6,
            )

        if ax_id == 2:
            label_text = (
                "Our approach - CI moment tensor\n"
                + r"$\Delta \vec{x} = $"
                + f"{dist:0.1f}km"
            )
            ax.set_xticks([-119, -118, -117])
            ax.set_xticklabels(
                [r"119$^\circ$W", r"118$^\circ$W", r"117$^\circ$W"], fontsize=6,
            )
            # cb_ax = fig.add_axes([0.915 + 0.015, 0.3, 0.015, 0.4])
            x0, y0, w, h = ax.get_position().bounds
            cb_ax = fig.add_axes([x0 + w + 0.05 * w, y0, 0.05 * w, h])
            cbar = plt.colorbar(pcm, cax=cb_ax)  # , extend="min")
            # cbar.add_lines(cnt)
            cbar.set_ticks([-1, 0, 1])
            cbar.set_label(label="Normalised Beampower", fontsize=6)
            cbar.ax.tick_params(labelsize=6)

        t = ax.text(
            0.05,
            0.95,
            label_text,
            ha="left",
            va="top",
            fontsize=6,
            # fontweight="bold",
            transform=ax.transAxes,
            zorder=100,
        )
        t.set_bbox(dict(facecolor="w", alpha=0.75, edgecolor="None"))

        lblbls = ["a)", "b)", "c)"]
        ax.text(
            0.00,
            1.01,
            lblbls[ax_id],
            ha="left",
            va="bottom",
            fontsize=10,
            # fontweight="bold",
            transform=ax.transAxes,
            zorder=100,
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
        # pcm = ax.contourf(xx, yy, beampowers.T, levels=np.linspace(-1, 1, 75), vmin=vmin, vmax=vmax, cmap=roma_map) # , norm=LogNorm())
        ax.set_xticks([-75, 0, 75])
        ax.set_xticklabels([-75, 0, 75], fontsize=6)
        ax.set_yticks([-75, 0, 75])
        ax.set_yticklabels([-75, 0, 75], fontsize=6)
        if ax_id == 0:
            if LHS:
                ax.set_ylabel("Distance [km]", fontsize=6)
                ax.set_xticklabels([])
                ax.set_yticklabels([-75, 0, 75], fontsize=6)
            else:
                ax.set_xticklabels([])
                ax.set_yticklabels([-75, 0, 75], fontsize=6)
        elif ax_id == 1:
            ax.set_xticklabels([])
            ax.set_yticklabels([])
        elif ax_id == 2:
            if LHS:
                ax.set_yticklabels([-75, 0, 75], fontsize=6)
                ax.set_ylabel("Distance [km]", fontsize=6)
                ax.set_xticklabels([-75, 0, 75], fontsize=6)
                # ax.set_xlabel("Distance [km]", fontsize=6)
            else:
                ax.set_yticklabels([-75, 0, 75], fontsize=6)
                ax.set_xticklabels([-75, 0, 75], fontsize=6)
                # ax.set_xlabel("Distance [km]", fontsize=6)
        elif ax_id == 3:
            ax.set_yticklabels([])
            # ax.set_xlabel("Distance [km]", fontsize=6)
            ax.set_xticklabels([-75, 0, 75], fontsize=6)
        if ax_id is None:
            ax.set_xlabel("Distance [km]", fontsize=6)
            if LHS:
                ax.set_ylabel("Distance [km]", fontsize=6)
        if eiwt_vel:
            ax.text(
                0.05,
                0.95,
                f"analytical GF, v={eiwt_vel}km/s",
                ha="left",
                va="top",
                fontsize=6,
                fontweight="bold",
                transform=ax.transAxes,
            )
        else:
            ax.text(
                0.05,
                0.95,
                f"Numerical GF",
                ha="left",
                va="top",
                fontsize=6,
                fontweight="bold",
                transform=ax.transAxes,
            )
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
    # if settings["do_synth"]:

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
            zorder=10,
        )
    elif settings["geometry_type"] == "cartesian":
        ax.scatter(
            lon_max,
            lat_max,
            facecolors="k",
            edgecolors="w",
            linewidth=0.5,
            s=50,
            marker="o",
            label="Beampower Peak",
        )

    if source_loc is not None:
        if settings["geometry_type"] == "geographic":
            ax.scatter(
                source_loc[0],
                source_loc[1],
                edgecolors="k",
                c="w",
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
    "/data/cen/u254/sven/global_mfp/projects/173_redo_chino_hills_denser/1636106117_settings.yml",
    "r",
) as stream:
    settings = yaml.safe_load(stream)
# fig, axs = plt.subplots(2, 2)
fig, axs = plt.subplots(
    1, 3, subplot_kw={"projection": ccrs.PlateCarree()}, figsize=(6, 1.875)
)
fig.subplots_adjust(left=0.06, right=0.9, top=0.95, wspace=0)

bp_fn = "/data/cen/u254/sven/global_mfp/projects/173_redo_chino_hills_denser_GF_15.5km/out/out_1217356935.0_3880_[0.1, 0.2]_0_0_False.npy"
# bp_fn = "/data/cen/u254/sven/global_mfp/projects/166_chino_GF_expl_MT_depth18km/out/out_1217356935.0_3880_[0.1, 0.5]_0_0_False.npy"
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
    ax=axs[1],
    vmin=-1,
    vmax=1,
    plot_station_locations=True,
    ax_id=1,
    LHS=True,
    fig=fig,
    # title="GF explosion",
)

# GF
bp_fn = "/data/cen/u254/sven/global_mfp/projects/173_redo_chino_hills_denser_GF_15.5km_trueGF/out/out_1217356935.0_3880_[0.1, 0.2]_0_0_False.npy"
# bp_fn = "/data/cen/u254/sven/global_mfp/projects/166_chino_GF_correct_MT_depth18km/out/out_1217356935.0_3880_[0.1, 0.5]_0_0_False.npy"
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
    ax=axs[2],
    ax_id=2,
    vmin=-1,
    vmax=1,
    plot_station_locations=True,
    LHS=True,
    fig=fig,
    # title="GF correct MT",
)

# vconst
bp_fn = "/data/cen/u254/sven/global_mfp/projects/173_redo_chino_hills_denser/out/out_1217356935.0_3880_[0.1, 0.2]_0_0_False.npy"
# bp_fn = "/data/cen/u254/sven/global_mfp/projects/166_chino_vconst3.0_depth18km/out/out_1217356935.0_3880_[0.1, 0.5]_0_0_False.npy"
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
    ax=axs[0],
    ax_id=0,
    vmin=-1,
    vmax=1,
    plot_station_locations=True,
    LHS=True,
    fig=fig,
    # title="v=3 km/s",
)

fig.savefig("/data/cen/u254/sven/global_mfp/figure_scripts/png/fig5.png", dpi=300)
