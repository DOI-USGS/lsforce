import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import xarray as xr
import cartopy.crs as ccrs
import warnings

KM_PER_M = 1 / 1000  # [km/m]


class LSTrajectory:
    """Class for force inversion derived trajectories.

    Attributes:
        force:
        mass_requested:
        target_length:
        duration:
        detrend_velocity:
    """

    def __init__(
        self,
        force,
        mass=None,
        target_length=None,
        duration=None,
        detrend_velocity=None,
    ):
        """
        Args:
            force (LSForce): Completed force inversion
            mass_requested:
            target_length:
            duration:
            detrend_velocity:
        """

        if not force.inversion_complete:
            raise ValueError('Cannot compute trajectory if inversion has not been run!')

        self.force = force
        self.mass_requested = mass
        self.target_length = target_length
        self.duration = duration
        self.detrend_velocity = detrend_velocity

        compute_kwargs = dict(
            mass=self.mass_requested,
            target_length=self.target_length,
            duration=self.duration,
            detrend_velocity=detrend_velocity,
            plot_jackknife=True,  # REMOVE!
        )

        self._compute_trajectory(**compute_kwargs)
        self._compute_trajectory(elevation_profile=True, **compute_kwargs)

    def _compute_trajectory(
        self,
        mass=None,
        target_length=None,
        duration=None,
        elevation_profile=False,
        plot_jackknife=False,
        image=None,
        dem=None,
        reference_point=None,
        detrend_velocity=None,
    ):
        """
        Integrate force time series to velocity and then displacement. Either
        provide a mass or a target horizontal runout length. If a length is
        provided, the code will find the mass that achieves this length. Calls
        _trajectory_automass().

        Args:
            mass: Landslide mass [kg]
            target_length: Horizontal runout length from groundtruth [m]
            duration: Clip time series to go from 0-duration [s]
            elevation_profile: If True, plot vertical displacement versus
                               horizontal runout distance (H vs. L) instead of
                               a map view
            plot_jackknife: Toggle plotting jackknifed displacements as well
            image: An xarray.DataArray with coordinates defined in km with the
                   origin (0, 0) being the start location of the compute_trajectory
            dem: A UTM-projected DEM GeoTIFF to slice thru for elevation
                 profile plot
            reference_point (int/float or list): Plot a dot on compute_trajectory, and
                                                 line on colorbar, at this
                                                 specified time(s) for
                                                 reference (default: None, for
                                                 no markings)
            detrend_velocity: If provided, force velocity to linearly go to
                              zero at this time [s]. If None, don't detrend

        Returns:
            The output figure
        """

        # Convert reference points to numpy array
        reference_points = np.atleast_1d(reference_point)

        # For the full inversion (all channels) result
        (
            self.Za,
            self.Ea,
            self.Na,
            self.Zvel,
            self.Evel,
            self.Nvel,
            self.Zdisp,
            self.Edisp,
            self.Ndisp,
            self.mass_actual,
            self.traj_tvec,
        ) = self._trajectory_automass(
            self.force.Zforce,
            self.force.Eforce,
            self.force.Nforce,
            mass=mass,
            target_length=target_length,
            duration=duration,
            detrend=detrend_velocity,
        )
        self.Hdist = _calculate_Hdist(self.Edisp, self.Ndisp)

        fig, ax = plt.subplots()

        # Converting to km below
        if elevation_profile:
            x = self.Hdist / 1000
            y = self.Zdisp / 1000
        else:
            x = self.Edisp / 1000
            y = self.Ndisp / 1000
        sc = ax.scatter(x, y, c=self.traj_tvec, cmap='rainbow', zorder=100)

        if elevation_profile:
            ax.set_xlabel('Horizontal distance (km)')
            ax.set_ylabel('Vertical distance (km)')
        else:
            ax.set_xlabel('East distance (km)')
            ax.set_ylabel('North distance (km)')

        t0 = self.force.st[0].stats.starttime
        if self.force.zero_time:
            t0 += self.force.zero_time
        cbar = plt.colorbar(
            sc, label='Time (s) from {}'.format(t0.strftime('%Y-%m-%d %H:%M:%S'))
        )

        # Plot reference points, if any
        if np.any(reference_points):
            cmap = cm.get_cmap('Greys_r', reference_points.size)
            for i, time in enumerate(reference_points):
                try:
                    ref_pt_ind = np.where(self.traj_tvec == time)[0][0]
                except IndexError:
                    raise  # No point corresponding to requested reference time
                ax.scatter(x[ref_pt_ind], y[ref_pt_ind], color=cmap(i), zorder=150)
                cbar.ax.plot(
                    [self.traj_tvec.min(), self.traj_tvec.max()],
                    [time, time],
                    color=cmap(i),
                    linewidth=2,
                )

        title = f'mass = {self.mass_actual:,} kg\nrunout length = {self.Hdist[-1] / 1000:.2f} km'
        if target_length:
            title += f'\n(target length = {target_length:g} km)'
        ax.set_title(title)

        # Plot jackknife trajectories as well if desired
        if plot_jackknife:
            self.jackknife = dict(num_iter=self.force.jackknife['num_iter'])
            self.jackknife['Zdisp_all'] = []
            self.jackknife['Edisp_all'] = []
            self.jackknife['Ndisp_all'] = []
            self.jackknife['Hdist_all'] = []
            for i in range(self.force.jackknife['num_iter']):
                *_, Zdisp_i, Edisp_i, Ndisp_i, _, _ = self._trajectory_automass(
                    self.force.jackknife['Zforce_all'][i],
                    self.force.jackknife['Eforce_all'][i],
                    self.force.jackknife['Nforce_all'][i],
                    mass=mass,
                    target_length=target_length,
                    duration=duration,
                    detrend=detrend_velocity,
                )

                # Store jackknifed trajectories
                self.jackknife['Zdisp_all'].append(Zdisp_i)
                self.jackknife['Edisp_all'].append(Edisp_i)
                self.jackknife['Ndisp_all'].append(Ndisp_i)

                # Converting to km below
                if elevation_profile:
                    Hdist_i = _calculate_Hdist(Edisp_i, Ndisp_i)
                    self.jackknife['Hdist_all'].append(Hdist_i)  # Store
                    x = Hdist_i / 1000
                    y = Zdisp_i / 1000
                else:
                    x = Edisp_i / 1000
                    y = Ndisp_i / 1000
                ax.scatter(x, y, c=self.traj_tvec, cmap='rainbow', alpha=0.02)

        ax.axis('equal')

        if (image is not None) and (not elevation_profile):
            xlim = ax.get_xlim()
            ylim = ax.get_ylim()
            image.plot.imshow(
                ax=ax, cmap='Greys_r', add_colorbar=False, add_labels=False, zorder=-10
            )
            ax.axis('equal')
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)

        if (dem is not None) and elevation_profile:
            distance, drop = self._slice_dem(dem)
            ax.plot(distance / 1000, drop / 1000, color='black', zorder=100)

        plt.tight_layout()
        plt.show()

        return fig

    def _integrate_acceleration(
        self, Zforce, Eforce, Nforce, mass, startidx, endidx, detrend=None
    ):

        traj_tvec = self.force.tvec[startidx : endidx + 1]

        dx = 1.0 / self.force.force_sampling_rate
        Za = -Zforce.copy()[startidx : endidx + 1] / mass
        Ea = -Eforce.copy()[startidx : endidx + 1] / mass
        Na = -Nforce.copy()[startidx : endidx + 1] / mass
        Zvel = np.cumsum(Za) * dx
        Evel = np.cumsum(Ea) * dx
        Nvel = np.cumsum(Na) * dx

        # Detrend is either None (no detrending) or a time where velo should
        # be fully tapered to zero
        if detrend:
            zeroidx = np.where(traj_tvec == detrend)[0][
                0
            ]  # Index corresponding to time where velo should be zero
            for comp in [Zvel, Evel, Nvel]:
                trend = np.linspace(0, comp[zeroidx], len(traj_tvec[:zeroidx]))
                comp[:zeroidx] -= trend
                comp[zeroidx:] = np.zeros(len(comp[zeroidx:]))

        Zdisp = np.cumsum(Zvel) * dx
        Edisp = np.cumsum(Evel) * dx
        Ndisp = np.cumsum(Nvel) * dx

        return Za, Ea, Na, Zvel, Evel, Nvel, Zdisp, Edisp, Ndisp, traj_tvec

    def _trajectory_automass(
        self,
        Zforce,
        Eforce,
        Nforce,
        mass=None,
        target_length=None,
        duration=None,
        detrend=None,
    ):
        """
        Calls _integrate_acceleration().
        """

        # Check args
        if mass and target_length:
            raise ValueError('Cannot specify both mass and target length!')
        if not mass and not target_length:
            raise ValueError('You must specify either mass or target length!')

        startidx = np.where(self.force.tvec == 0)[0][0]  # Always start at t = 0
        # Clip time series to `duration` [s] if desired
        if duration:
            endidx = np.where(self.force.tvec == duration)[0][0]
        else:
            endidx = len(self.force.tvec)

        # Either use the mass that was provided, or calculate one
        if target_length:

            MASS_INC = int(1e7)  # [kg] Smaller increment is slower but more precise

            # Initialize with end-members
            mass = 0  # [kg]
            current_length = np.inf  # [km]

            while current_length > target_length:

                mass += MASS_INC  # Increase the mass

                # Calculate the runout length [km] based on this mass
                *_, Edisp, Ndisp, _ = self._integrate_acceleration(
                    Zforce, Eforce, Nforce, mass, startidx, endidx, detrend
                )
                current_length = _calculate_Hdist(Edisp, Ndisp)[-1] / 1000  # [km]
        else:
            mass = int(mass)

        # Calculate compute_trajectory based on mass assigned above
        (
            Za,
            Ea,
            Na,
            Zvel,
            Evel,
            Nvel,
            Zdisp,
            Edisp,
            Ndisp,
            traj_tvec,
        ) = self._integrate_acceleration(
            Zforce, Eforce, Nforce, mass, startidx, endidx, detrend
        )

        return Za, Ea, Na, Zvel, Evel, Nvel, Zdisp, Edisp, Ndisp, mass, traj_tvec

    def _slice_dem(self, dem_file, interp_spacing=0.1):
        """
        Slice through an input DEM along the compute_trajectory path.

        Args:
            dem_file: DEM GeoTIFF to slice (must be UTM-projected!)
            interp_spacing: [m] Density of interpolation points

        Returns:
            horizontal_distance: Distance along path given by EDisp, NDisp
            elevation: Elevation along horizontal_distance
        """

        dem = xr.open_rasterio(dem_file).squeeze()
        # Set no data values to NaN
        dem = dem.where(dem != dem.nodatavals)

        # Define interpolation points in UTM space
        crs = ccrs.epsg(int(dem.crs.split(':')[-1]))
        if crs.proj4_params['proj'] != 'utm':
            raise ValueError('Input DEM must have a UTM projection!')
        loc_utm = crs.transform_point(self.force.lon, self.force.lat, ccrs.Geodetic())
        points = [
            [x + loc_utm[0], y + loc_utm[1]] for x, y in zip(self.Edisp, self.Ndisp)
        ]

        # Densify the coarse points
        path_x, path_y = [], []
        x_prev, y_prev = points[0]
        for pt in points[1:]:
            x, y = pt

            seg_length = np.linalg.norm([y - y_prev, x - x_prev])
            # Choose a number of pts that gives approximately interp_spacing
            n = int(seg_length / interp_spacing)

            # Append densified path
            path_x = np.hstack([path_x, np.linspace(x_prev, x, n)])
            path_y = np.hstack([path_y, np.linspace(y_prev, y, n)])

            x_prev, y_prev = x, y

        # Actually interpolate!
        profile = dem.interp(
            x=xr.DataArray(path_x), y=xr.DataArray(path_y), method='linear'
        )

        # Find horizontal distance vector (works for curvy paths!)
        horiz_dist = np.hstack(
            [
                0,
                np.cumsum(
                    np.linalg.norm([np.diff(profile.x), np.diff(profile.y)], axis=0)
                ),
            ]
        )

        # Check that interp_spacing wasn't too coarse by matching path lengths
        if not np.isclose(horiz_dist[-1], self.Hdist[-1]):
            raise ValueError('interp_spacing was too coarse. Try decreasing.')

        warnings.warn('Assuming DEM vertical unit is meters!')

        profile.data -= profile.data[0]  # Start at 0 and go negative

        return horiz_dist, profile.data


def _calculate_Hdist(Edisp, Ndisp):
    """
    Calculate horizontal distance vector (Hdist) from east and north
    displacement vectors. This is the horizontal distance "along the avalanche
    path" as a function of time. Hdist[-1] is L, the horizontal runout distance
    (which is shorter than the 3-D runout distance).

    Args:
        Edisp: Eastward displacement vector as a function of time [m]
        Ndisp: Northward displacement vector as a function of time [m]

    Returns:
        Horizontal distance as a function of time [m]
    """

    dx = np.diff(Edisp)
    dy = np.diff(Ndisp)

    return np.hstack([0, np.cumsum(np.linalg.norm([dx, dy], axis=0))])
