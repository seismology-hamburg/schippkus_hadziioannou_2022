import datetime
import os
from itertools import product

import cartopy.crs as ccrs
import numpy as np
import obspy
import pylab as plt
import xarray
import yaml
from astropy.convolution import Gaussian2DKernel, convolve
from cmcrameri import cm
from global_land_mask import is_land
from tqdm import tqdm


def get_plotting_info_for_project(settings):
    import obspy

    def get_station_locs_from_staxml(st, settings):
        import obspy

        station_locations = []
        for tr in st:
            inv = obspy.read_inventory(
                f"/data/cen/u254/sven/global_mfp/data/stations/{tr.stats.network}.{tr.stats.station}.xml",
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
    # st = obspy.read("./data/deconv/CI/LHZ_chino.mseed")
    # st = obspy.read("./data/deconv/northsea_2019_02/*.mseed")
    # st.merge(fill_value=0)
    # print(st)
    # # station_locations = get_station_locs_from_staxml(st, settings)
    station_locations = np.load(
        "/data/cen/u254/sven/global_mfp/figure_scripts/fig6_stations.npy"
    )
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
    triangle_size=4,
    zoom=None,
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
        ax.coastlines(resolution="10m", linewidth=0.5, color="k")

    if plot_station_locations:
        if settings["geometry_type"] == "geographic":
            ax.scatter(
                station_locations[:, 0],
                station_locations[:, 1],
                c="k",
                s=triangle_size,
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
            cmap=cm.vik_r,
        )
        if eiwt_vel:
            lblx, lbly = 0.05, 0.95
            if "km/s" in eiwt_vel:
                lblx, lbly = 0.07, 0.93
            t = ax.text(
                lblx,
                lbly,
                f"{eiwt_vel}",
                ha="left",
                va="top",
                fontsize=8,
                transform=ax.transAxes,
            )
            t.set_bbox(dict(facecolor="w", alpha=0.75, edgecolor="None"))
        # ax.set_xlim(
        #     settings["grid_limits_lon"][0] + 0.75, settings["grid_limits_lon"][1] - 1
        # )
        # ax.set_ylim(settings["grid_limits_lat"][0], settings["grid_limits_lat"][1])
        # ax.set_xlim(source_loc[0] - 2, source_loc[0] + 2)
        # ax.set_ylim(source_loc[1] - 1.5, source_loc[1] + 1.5)
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
            t = ax.text(
                0.05,
                0.95,
                f"{eiwt_vel}",
                ha="left",
                va="top",
                fontsize=6,
                transform=ax.transAxes,
            )
            t.set_bbox(dict(facecolor="w", alpha=0.75, edgecolor="None"))
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
    max_indices = np.unravel_index(np.argmax(beampowers, axis=None), beampowers.shape)
    lon_max = lons[max_indices[0]]
    lat_max = lats[max_indices[1]]
    # if settings["geometry_type"] == "geographic":
    #     ax.scatter(
    #         lon_max,
    #         lat_max,
    #         facecolors="k",
    #         edgecolors="w",
    #         linewidth=0.5,
    #         s=50,
    #         marker="o",
    #         label="Beampower Peak",
    #         transform=trans,
    #     )
    # elif settings["geometry_type"] == "cartesian":
    #     ax.scatter(
    #         lon_max,
    #         lat_max,
    #         facecolors="k",
    #         edgecolors="w",
    #         linewidth=0.5,
    #         s=50,
    #         marker="o",
    #         label="Beampower Peak",
    #     )

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
            print(
                obspy.geodetics.gps2dist_azimuth(
                    lat1=source_loc[1], lon1=source_loc[0], lat2=lat_max, lon2=lon_max
                )[0]
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


# project_id = '118_northsea_one_day'
# project_id = '124_northsea_rerun_oneday_fixednorm'
# project_id = '128_atlantic_oneday'
# project_id = '128_northsea_oneday_secondary'
# project_id = "131_northsea_2019_02"
project_id_1 = "161_2019_02_GF_one_week_1min"
project_id_2 = "161_2019_02_vconst2.6_one_week_1min"
# project_id = '129_northsea_2014_02_outlier_in_logic'
# project_id = '134_northsea_2019_02_vconst3.5'
# project_id = '149_northsea_multifreq'
# project_id = '149_northsea_multifreq_KESW_elim'
month_used = "2019-02"

project_ids = [
    # "161_2019_02_GF_one_week_1min_mediumgrid_rerun",
    # "169_redo_2019_02_medium_one_week_GF",
    # "169_redo_2019_02_medium_one_week_vconst2.6",
    # "169_redo_2019_02_medium_one_week_vconst2.8",
    # "169_redo_2019_02_medium_one_week_vconst3.0",
    # "171_redo_2019_02_medium_one_week_GF",
    # "173_redo_2019_02_one_week_v3.2",
    "174_redo_2019_02_GF",
    "174_redo_2019_02_v3.2",
    # "169_redo_2019_02_medium_one_week_vconst3.2",
    # "161_2019_02_GF_one_week_1min_largegrid_oldcsdm",
    # "161_2019_02_vconst2.6_one_week_1min",
    # "161_2019_02_vconst2.8_one_week_1min_qual",
    # "161_2019_02_vconst3.0_one_week_1min_largegrid",
    # "161_2019_02_vconst3.2_one_week_1min_qual",
]

wdirs = [
    f"/data/cen/u254/sven/global_mfp/projects/{project_id}/"
    for project_id in project_ids
]


# with open(
#     f"/data/cen/u254/sven/global_mfp/projects/161_2019_02_GF_one_week_1min/1633012866.0659773_settings.yml",
#     "r",
# ) as stream:
#     settings = yaml.safe_load(stream)

measure = "p2l"
measure = "hs"

# nc_fn = '../data/wavewatch/WW3_GLOB_30M_201902.nc'
# nc_fn = '../data/wavewatch/WW3-ATNE-10M_201902_p2l.nc'
nc_fn = "/data/cen/u254/sven/global_mfp/data/wavewatch/WW3-GLOB-30M_201902_p2l.nc"
nc_fn = "/data/cen/u254/sven/global_mfp/data/wavewatch/WW3-GLOB-30M_201902_hs.nc"

# filterbands = ["[0.065, 0.075]", "[0.13, 0.15]"]
filterbands = ["[0.13, 0.15]"]
# filterbands = ["_"]


settings_list = []
for wdir in wdirs:
    # get latest settings file
    highest_timestamp = 0
    for settings_file in os.listdir(f"{wdir}/"):
        if ".yml" in settings_file:
            timestamp = obspy.UTCDateTime(float(settings_file.split("_")[0]))
            if timestamp > highest_timestamp:
                highest_timestamp = timestamp
    # if "v3.2" in wdir:
    with open(f"{wdir}/{int(highest_timestamp.timestamp)}_settings.yml", "r") as stream:
        settings_list.append(yaml.safe_load(stream))
    # else:
    #     with open(f"{wdir}/{highest_timestamp.timestamp}_settings.yml", "r") as stream:
    #         settings_list.append(yaml.safe_load(stream))

print("read ww3 data")
# nc_fn = '../data/wavewatch/WW3_GLOB_30M_201402_hs.nc'
# nc_fn = '../data/wavewatch/WW3_GLOB_30M_201405_p2l.nc'
dataset = xarray.open_dataset(nc_fn, decode_times=True)

print(dataset)


def time_is_in_timewindow(current_time, filetime):
    # 3 hr steps in WW3 data
    timestamp_current = current_time.data.astype(float) / 1e9
    timestamp_filetime = float(filetime)
    if (
        timestamp_current - 1.5 * 3600
        <= timestamp_filetime
        <= timestamp_current + 1.5 * 3600
    ):
        return True
    return False


# gps_are_land = is_land(grid_points[:, 1], grid_points[:, 0])

# mean_bps_pertime = []
# ww_per_time = []
bps_per_time_and_project = []
xx_yy_SWH = []
swhs = []
used_times = []
for current_time, filterband in tqdm(
    product(dataset.sel(time=month_used).time, filterbands)
):
    # select times to plot here
    if current_time.data.astype(float) / 1e9 in [
        1549065600.0,
        1549173600.0,
        1549465200.0,
    ]:
        print(current_time)
    else:
        continue
    # current_time.data.astype(int)/1E9 > obspy.UTCDateTime('2019-02-28T12:00:00.0Z').timestamp
    used_times.append(current_time)

    bps_per_project = []
    for wdir in wdirs:
        beampowers = []
        for filename in os.listdir(f"{wdir}/out/"):
            if "out_" in filename:
                filetime = obspy.UTCDateTime(float(filename.split("_")[1]))
                if time_is_in_timewindow(current_time, filetime):
                    bp = np.load(f"{wdir}/out/{filename}")
                    bp /= np.nanmax(np.abs(bp))
                    # bp_conv = convolve(bp, kernel)
                    beampowers.append(bp)

        if len(beampowers) == 0:
            print(f"no data for {current_time.data.astype(float)/1E9}")
            continue

        beampowers = np.array(beampowers)
        mean_bp = np.mean(beampowers, axis=0)

        bps_per_project.append(mean_bp)
    bps_per_time_and_project.append(bps_per_project)

    # print(dataset)
    if measure == "p2l":
        # print(dataset.p2l.sel(time='2014-05-23').shape)
        # print(dataset.p2l.sel(time='2014-05-23').mean(axis=0).sel(f=0.125, method='nearest'))
        # ds_day = dataset.p2l.sel(time='2014-05-23')
        # mean_in_time = ds_day.sel(f=0.125, method='nearest').mean(axis=0)
        # print(mean_in_time.shape)
        # hs_mean_in_time = dataset.p2l.sel(time='2014-05-23')[:,:,:].mean(axis=0)
        if filterband == "[0.065, 0.075]":
            f = 0.07
        elif filterband == "[0.13, 0.15]":
            f = 0.14
        wavewatch_data = dataset.p2l.sel(f=f, method="nearest")
        wavewatch_data = wavewatch_data.sel(time=current_time)
        # mean_in_time = ds_day.sel(f=f, method='nearest').mean(axis=0)
        xx, yy = np.meshgrid(dataset.p2l.longitude, dataset.p2l.latitude)
    if measure == "hs":
        # mean_in_time = dataset.hs.sel(time='2019-09-24')[:,:,:].mean(axis=0)
        wavewatch_data = dataset.hs.sel(time=current_time)
        xx, yy = np.meshgrid(dataset.hs.longitude, dataset.hs.latitude)

    swhs.append(wavewatch_data)
    xx_yy_SWH.append([xx, yy])

subplot_kw = {"projection": ccrs.PlateCarree()}
fig, axs = plt.subplots(3, 3, subplot_kw=subplot_kw, figsize=(8, 7))
fig.subplots_adjust(
    left=0.01, right=0.99, wspace=0, hspace=0.075, bottom=0.1, top=0.975
)
eu_extent = [
    settings_list[0]["grid_limits_lon"][0],
    settings_list[0]["grid_limits_lon"][1],
    settings_list[0]["grid_limits_lat"][0],
    settings_list[0]["grid_limits_lat"][1] - 0.5,
]

from shapely.geometry import box

zoom_lims = [2, 9, eu_extent[2] + 0.5, eu_extent[2] + 4.5]
zoom_box = box(zoom_lims[0], zoom_lims[2], zoom_lims[1], zoom_lims[3])

idx_map = {0: 1, 1: 2, 2: 0}
figlabel_map = {
    0: ["a)", "b)", "c)"],
    1: ["d)", "e)", "f)"],
    2: ["g)", "h)", "i)"],
}

outer_grid = fig.add_gridspec(3, 3, wspace=0)
for idx, (bps_per_proj, time, [xx, yy], swh_map) in enumerate(
    zip(bps_per_time_and_project, used_times, xx_yy_SWH, swhs)
):
    # map idx to ax_idx:
    axidx = idx_map[idx]

    if len(bps_per_proj) != len(project_ids):
        print(f"skipping {time}")
        continue
    print("read waveform data")
    print(len(bps_per_proj))

    # left plots (broadband)
    # e-iwt - subplots for vconst=2.5, 2.6, 2.7, 2.8
    inner_grid_left = outer_grid[axidx, 0].subgridspec(2, 2, wspace=0, hspace=0)
    # axs_for_each_vconst = inner_grid_left.subplots(projection=ccrs.PlateCarree())
    gridspec_idxs = [[0, 0], [0, 1], [1, 0], [1, 1]]
    eiwts = [2.6, 2.8, 3.0, 3.2]
    eiwts = [3.0]
    quad_axs = []
    for settings, bps, gridspec_idx, eiwt_vel in zip(
        settings_list[1:], bps_per_proj[1:], gridspec_idxs, eiwts
    ):
        # curr_ax = fig.add_subplot(
        #     inner_grid_left[gridspec_idx[0], gridspec_idx[1]],
        #     projection=ccrs.PlateCarree(),
        # )
        curr_ax = axs[axidx, 0]
        quad_axs.append(curr_ax)
        (
            grid_lon_coords,
            grid_lat_coords,
            station_locations,
            source_loc,
        ) = get_plotting_info_for_project(settings)

        try:
            pcm_bp = plot_beampowers_on_map(
                lons=grid_lon_coords,
                lats=grid_lat_coords,
                beampowers=bps,
                settings=settings,
                outfile=None,
                source_loc=None,
                # title="v=2.6 km/s",
                fig=fig,
                ax=curr_ax,
                vmin=-1,
                vmax=1,
                station_locations=station_locations,
                plot_station_locations=True,
                triangle_size=2,
                # eiwt_vel=f"v={eiwt_vel}km/s",
                eiwt_vel=f"Standard approach",
            )
        except TypeError:
            print(settings["project_id"])
            print(grid_lon_coords.shape)
            print(bps.shape)
        curr_ax.set_extent(eu_extent)

        if gridspec_idx == [0, 0]:
            curr_ax.text(
                0.00,
                1.01,
                figlabel_map[axidx][0],
                ha="left",
                va="bottom",
                fontsize=12,
                # fontweight="bold",
                transform=curr_ax.transAxes,
                zorder=100,
            )

        if None:
            # if idx == 1 and gridspec_idx == [0, 0]:
            zoom_ax_2 = fig.add_axes(
                [0.025, 0.1, 0.15, 0.15], projection=ccrs.PlateCarree()
            )
            pcm = plot_beampowers_on_map(
                lons=grid_lon_coords,
                lats=grid_lat_coords,
                beampowers=bps,
                settings=settings_list[0],
                outfile=None,
                source_loc=None,
                # title="GF",
                fig=fig,
                ax=zoom_ax_2,
                vmin=-1,
                vmax=1,
                plot_station_locations=False,
                eiwt_vel=None,
            )
            zoom_ax_2.set_extent(zoom_lims)
            zoom_ax_2.add_geometries(
                [zoom_box], crs=ccrs.PlateCarree(), facecolor="None", edgecolor="k"
            )

    (
        grid_lon_coords,
        grid_lat_coords,
        station_locations,
        source_loc,
    ) = get_plotting_info_for_project(settings_list[0])
    # print(settings_list[0])
    pcm = plot_beampowers_on_map(
        lons=grid_lon_coords,
        lats=grid_lat_coords,
        beampowers=bps_per_proj[0],
        settings=settings_list[0],
        outfile=None,
        source_loc=None,
        # title="GF",
        fig=fig,
        ax=axs[axidx, 1],
        vmin=-1,
        vmax=1,
        station_locations=station_locations,
        plot_station_locations=True,
        triangle_size=2,
        eiwt_vel="Our approach",
    )

    if None:
        # if idx == 1:
        zoom_ax_1 = fig.add_axes([0.35, 0.1, 0.15, 0.15], projection=ccrs.PlateCarree())
        pcm = plot_beampowers_on_map(
            lons=grid_lon_coords,
            lats=grid_lat_coords,
            beampowers=bps_per_proj[0],
            settings=settings_list[0],
            outfile=None,
            source_loc=None,
            # title="GF",
            fig=fig,
            ax=zoom_ax_1,
            vmin=-1,
            vmax=1,
            plot_station_locations=False,
            eiwt_vel=None,
        )
        zoom_ax_1.set_extent(zoom_lims)
        zoom_ax_shw = fig.add_axes(
            [0.675, 0.1, 0.15, 0.15], projection=ccrs.PlateCarree()
        )
        pcm = zoom_ax_shw.pcolormesh(
            xx,
            yy,
            swh_map,
            edgecolors="face",
            transform=ccrs.PlateCarree(),
            vmin=-1,
            vmax=8,
            cmap=cm.lapaz,
        )

        zoom_ax_shw.coastlines(resolution="10m", linewidth=0.5, color="k")
        zoom_ax_shw.set_extent(zoom_lims)

    axs[axidx, 1].set_extent(eu_extent)

    axs[axidx, 1].text(
        0.00,
        1.01,
        figlabel_map[axidx][1],
        ha="left",
        va="bottom",
        fontsize=12,
        # fontweight="bold",
        transform=axs[axidx, 1].transAxes,
        zorder=100,
    )
    # t = axs[axidx, 1].text(
    #     0.05,
    #     0.95,
    #     "numerical GF",
    #     ha="left",
    #     va="top",
    #     fontsize=8,
    #     # fontweight="bold",
    #     transform=axs[axidx, 1].transAxes,
    #     zorder=100,
    # )
    # t.set_bbox(dict(facecolor="w", alpha=0.75, edgecolor="None"))

    # # fix subaxes y-position
    # tl_x0, tl_y0, w, h = quad_axs[0].get_position().bounds
    # tr_x0, tr_y0, w, h = quad_axs[1].get_position().bounds
    # bl_x0, bl_y0, w, h = quad_axs[2].get_position().bounds
    # br_x0, br_y0, w, h = quad_axs[3].get_position().bounds
    # vertical_space = tl_y0 - (bl_y0 + h)
    # quad_axs[0].set_position(pos=[tl_x0, tl_y0 - (vertical_space / 2), w, h])
    # quad_axs[1].set_position(pos=[tr_x0, tr_y0 - (vertical_space / 2), w, h])
    # quad_axs[2].set_position(pos=[bl_x0, bl_y0 + (vertical_space / 2), w, h])
    # quad_axs[3].set_position(pos=[br_x0, br_y0 + (vertical_space / 2), w, h])

    ax = axs[axidx, 2]
    pcm = ax.pcolormesh(
        xx,
        yy,
        swh_map,
        edgecolors="face",
        transform=ccrs.PlateCarree(),
        vmin=-1,
        vmax=8,
        cmap=cm.lapaz,
    )

    ax.coastlines(resolution="10m", linewidth=0.5, color="k")
    ax.set_extent(eu_extent)

    axs[axidx, 2].text(
        0.00,
        1.01,
        figlabel_map[axidx][2],
        ha="left",
        va="bottom",
        fontsize=12,
        # fontweight="bold",
        transform=axs[axidx, 2].transAxes,
        zorder=100,
    )

    time_string = datetime.datetime.fromtimestamp(
        int(time.data.astype(int) / 1e9)
    ).strftime("%Y-%m-%d %Hh")
    t = axs[axidx, 2].text(
        0.05,
        0.95,
        f"Significant wave height",
        ha="left",
        va="top",
        fontsize=8,
        # fontweight="bold",
        transform=axs[axidx, 2].transAxes,
        zorder=100,
    )
    t.set_bbox(dict(facecolor="w", alpha=0.75, edgecolor="None"))
    axs[axidx, 1].text(
        0.5,
        1.00,
        f"{time_string}",
        ha="center",
        va="bottom",
        fontsize=12,
        # fontweight="bold",
        transform=axs[axidx, 1].transAxes,
        zorder=100,
    )

    # t = ax.text(
    #     0.05,
    #     0.95,
    #     "Significant Wave Height",
    #     ha="left",
    #     va="top",
    #     fontsize=8,
    #     # fontweight="bold",
    #     transform=ax.transAxes,
    #     zorder=100,
    # )
    # t.set_bbox(dict(facecolor="w", alpha=0.75, edgecolor="None"))

    # ax.set_title(
    #     datetime.datetime.fromtimestamp(
    #         int(current_time.data.astype(int) / 1e9)
    #     ).strftime("%Y-%m-%d %H:00:00")
    # )

if None:
    for _ax in axs[2, :]:
        _ax.add_geometries(
            [zoom_box],
            crs=ccrs.PlateCarree(),
            facecolor="None",
            edgecolor="k",
            lw=2,
            zorder=10,
        )

x0, y0, w, h = axs[2, 0].get_position().bounds
x0, _, w, _ = axs[2, 2].get_position().bounds
# cb_ax = fig.add_axes([x0 + w + 0.025 * w, y0, 0.025 * w, h])
cb_ax = fig.add_axes([x0 - 1.5 * w, y0 - 2 * 0.035 * h, w, 0.035 * h])
cb = plt.colorbar(pcm_bp, cax=cb_ax, orientation="horizontal", extend="both")
if measure == "hs":
    cb.set_label("Normalized Beampower")
if measure == "p2l":
    cb.set_label(r"P2L [Pa$^2$m$^2$s]")

x0, y0, w, h = axs[2, 2].get_position().bounds
# cb_ax = fig.add_axes([x0 + w + 0.025 * w, y0, 0.025 * w, h])
cb_ax = fig.add_axes([x0, y0 - 2 * 0.035 * h, w, 0.035 * h])
cb = plt.colorbar(pcm, cax=cb_ax, orientation="horizontal")
if measure == "hs":
    cb.set_label("Significant Wave Height [m]")
if measure == "p2l":
    cb.set_label(r"P2L [Pa$^2$m$^2$s]")

# print(datetime.datetime(current_time.data.astype(float)))
# datetime.fromtimestamp(current_time.data)

# axs[0, 0].set_title("Standard approach")
# axs[0, 1].set_title("Our approach")
# axs[0, 2].set_title("Significant Wave Height")

# ax.set_extent([-30, 20, 30, 70])
fig.savefig(
    f"/data/cen/u254/sven/global_mfp/figure_scripts/png/fig6_5.png", dpi=300,
)
plt.close(fig)
# animation
