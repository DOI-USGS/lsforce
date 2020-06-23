import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import xarray as xr
import cartopy.crs as ccrs
from obspy.core import AttribDict
import warnings

KM_PER_M = 1 / 1000  # [km/m]

MASS_INC = int(1e7)  # [kg] Smaller increment is slower but more precise


class LSTrajectory:
    """Class for force inversion derived trajectories.

    Attributes:
        force:
        mass_requested:
        target_length:
        duration:
        detrend_velocity:
        jackknife:
        acceleration:
        velocity:
        displacement:
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
            mass:
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
        self.jackknife = None

        compute_kwargs = dict(
            mass=self.mass_requested,
            target_length=self.target_length,
            duration=self.duration,
            detrend_velocity=detrend_velocity,
        )

        self._compute_trajectory(**compute_kwargs)

    def plot_trajectory(
        self,
        elevation_profile=False,
        plot_jackknife=False,
        image=None,
        dem=None,
        reference_point=None,
    ):
        """Plot trajectory results with context.

        Args:
            elevation_profile: If True, plot vertical displacement versus
                               horizontal runout distance (H vs. L) instead of
                               a map view
            plot_jackknife: Toggle plotting jackknifed displacements as well (if
                available)
            image: An xarray.DataArray with coordinates defined in km with the
                   origin (0, 0) being the start location of the compute_trajectory
            dem: A UTM-projected DEM GeoTIFF to slice thru for elevation
                 profile plot
            reference_point (int/float or list): Plot a dot on compute_trajectory, and
                                                 line on colorbar, at this
                                                 specified time(s) for
                                                 reference (default: None, for
                                                 no markings)

        Returns:
            The output figure
        """

        # Convert reference points to numpy array
        reference_points = np.atleast_1d(reference_point)

        fig, ax = plt.subplots()

        # Converting to km below
        if elevation_profile:
            x = self.horizontal_distance * KM_PER_M
            y = self.displacement.Z * KM_PER_M
        else:
            x = self.displacement.E * KM_PER_M
            y = self.displacement.N * KM_PER_M
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

        title = (
            f'mass = {self.mass_actual:,} kg\n'
            f'runout length = {self.horizontal_distance[-1] * KM_PER_M:.2f} km'
        )
        if self.target_length:
            title += f'\n(target length = {self.target_length:g} km)'
        ax.set_title(title)

        # Plot jackknife trajectories as well if desired
        if plot_jackknife:
            if self.jackknife:
                for i in range(self.force.jackknife.num_iter):
                    if elevation_profile:
                        x = self.jackknife.horizontal_distance[i] * KM_PER_M
                        y = self.jackknife.displacement.Z[i] * KM_PER_M
                    else:
                        x = self.jackknife.displacement.E[i] * KM_PER_M
                        y = self.jackknife.displacement.N[i] * KM_PER_M
                    ax.scatter(x, y, c=self.traj_tvec, cmap='rainbow', alpha=0.02)
            else:
                warnings.warn('No jackknife iterations to plot.')

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
            ax.plot(distance * KM_PER_M, drop * KM_PER_M, color='black', zorder=100)

        fig.tight_layout()
        fig.show()

        return fig

    def _compute_trajectory(
        self, mass=None, target_length=None, duration=None, detrend_velocity=None,
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
            detrend_velocity: If provided, force velocity to linearly go to
                              zero at this time [s]. If None, don't detrend
        """

        # For the full inversion (all channels) result
        (
            self.acceleration,
            self.velocity,
            self.displacement,
            self.mass_actual,
            self.traj_tvec,
        ) = self._trajectory_automass(
            self.force.Z,
            self.force.E,
            self.force.N,
            mass=mass,
            target_length=target_length,
            duration=duration,
            detrend=detrend_velocity,
        )
        self.horizontal_distance = _calculate_horizontal_distance(
            self.displacement.E, self.displacement.N
        )

        # Compute jackknife trajectories as well if the inversion was jackknifed
        if self.force.jackknife:
            self.jackknife = AttribDict(
                num_iter=self.force.jackknife.num_iter,
                displacement=AttribDict(Z=[], E=[], N=[]),
                horizontal_distance=[],
            )
            for i in range(self.jackknife.num_iter):
                *_, disp_i, _, _ = self._trajectory_automass(
                    self.force.jackknife.Z.all[i],
                    self.force.jackknife.E.all[i],
                    self.force.jackknife.N.all[i],
                    mass=mass,
                    target_length=target_length,
                    duration=duration,
                    detrend=detrend_velocity,
                )

                horiz_dist_i = _calculate_horizontal_distance(disp_i.E, disp_i.N)

                # Store jackknifed trajectories
                self.jackknife.displacement.Z.append(disp_i.Z)
                self.jackknife.displacement.E.append(disp_i.E)
                self.jackknife.displacement.N.append(disp_i.N)
                self.jackknife.horizontal_distance.append(horiz_dist_i)

    def _integrate_acceleration(
        self, z_force, e_force, n_force, mass, startidx, endidx, detrend=None
    ):

        traj_tvec = self.force.tvec[startidx : endidx + 1]

        dx = 1.0 / self.force.force_sampling_rate
        acceleration = AttribDict(
            Z=-z_force.copy()[startidx : endidx + 1] / mass,
            E=-e_force.copy()[startidx : endidx + 1] / mass,
            N=-n_force.copy()[startidx : endidx + 1] / mass,
        )
        velocity = AttribDict(
            Z=np.cumsum(acceleration.Z) * dx,
            E=np.cumsum(acceleration.E) * dx,
            N=np.cumsum(acceleration.N) * dx,
        )

        # Detrend is either None (no detrending) or a time where velocity should
        # be fully tapered to zero
        if detrend:
            # Index corresponding to time where velocity should be zero (closest entry)
            zeroidx = (np.abs(traj_tvec - detrend)).argmin()
            for comp in velocity.values():
                trend = np.linspace(0, comp[zeroidx], len(traj_tvec[:zeroidx]))
                comp[:zeroidx] -= trend
                comp[zeroidx:] = np.zeros(len(comp[zeroidx:]))

        displacement = AttribDict(
            Z=np.cumsum(velocity.Z) * dx,
            E=np.cumsum(velocity.E) * dx,
            N=np.cumsum(velocity.N) * dx,
        )

        return (
            acceleration,
            velocity,
            displacement,
            traj_tvec,
        )

    def _trajectory_automass(
        self,
        z_force,
        e_force,
        n_force,
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

        # Start as close to t = 0 as possible
        startidx = (np.abs(self.force.tvec - 0)).argmin()

        # Clip time series as close to `duration` [s] as possible, if desired
        if duration:
            endidx = (np.abs(self.force.tvec - duration)).argmin()
        else:
            endidx = len(self.force.tvec)

        # Either use the mass that was provided, or calculate one
        if target_length:

            # Initialize with end-members
            mass = 0  # [kg]
            current_length = np.inf  # [km]

            while current_length > target_length:

                mass += MASS_INC  # Increase the mass

                # Calculate the runout length [km] based on this mass
                *_, disp, _ = self._integrate_acceleration(
                    z_force, e_force, n_force, mass, startidx, endidx, detrend
                )
                current_length = (
                    _calculate_horizontal_distance(disp.E, disp.N)[-1] * KM_PER_M
                )  # [km]
        else:
            mass = int(mass)

        # Calculate compute_trajectory based on mass assigned above
        (
            acceleration,
            velocity,
            displacement,
            traj_tvec,
        ) = self._integrate_acceleration(
            z_force, e_force, n_force, mass, startidx, endidx, detrend
        )

        return (
            acceleration,
            velocity,
            displacement,
            mass,
            traj_tvec,
        )

    def _slice_dem(self, dem_file, interp_spacing=0.1):
        """
        Slice through an input DEM along the compute_trajectory path.

        Args:
            dem_file: DEM GeoTIFF to slice (must be UTM-projected!)
            interp_spacing: [m] Density of interpolation points

        Returns:
            horizontal_distance: Distance along path
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
            [x + loc_utm[0], y + loc_utm[1]]
            for x, y in zip(self.displacement.E, self.displacement.N)
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
        horizontal_distance = np.hstack(
            [
                0,
                np.cumsum(
                    np.linalg.norm([np.diff(profile.x), np.diff(profile.y)], axis=0)
                ),
            ]
        )

        # Check that interp_spacing wasn't too coarse by matching path lengths
        if not np.isclose(horizontal_distance[-1], self.horizontal_distance[-1]):
            raise ValueError('interp_spacing was too coarse. Try decreasing.')

        warnings.warn('Assuming DEM vertical unit is meters!')

        profile.data -= profile.data[0]  # Start at 0 and go negative

        return horizontal_distance, profile.data


def _calculate_horizontal_distance(east_displacement, north_displacement):
    """
    Calculate horizontal distance vector (horizontal_distance) from east and north
    displacement vectors. This is the horizontal distance "along the avalanche path" as
    a function of time. horizontal_distance[-1] is L, the horizontal runout distance
    (which is shorter than the 3-D runout distance).

    Args:
        east_displacement: Eastward displacement vector as a function of time [m]
        north_displacement: Northward displacement vector as a function of time [m]

    Returns:
        Horizontal distance as a function of time [m]
    """

    dx = np.diff(east_displacement)
    dy = np.diff(north_displacement)

    return np.hstack([0, np.cumsum(np.linalg.norm([dx, dy], axis=0))])
