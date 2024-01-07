import numpy as np


def my_target_lightcurve(data, ax=None, fig=None, index=None, zp=25,
                                lc_prop={}, bands=None, show_truth=True,
                                format_time=True, t0_format="mjd", 
                                phase_window=None, **kwargs):
        """ if index is None, a random index will be used. 
        if bands is None, the target's observed band will be used.
        """
        from matplotlib.colors import to_rgba
        from skysurvey.config import get_band_color
        
        if format_time:
            from astropy.time import Time
        if index is None:
            index = np.random.choice(data.obs_index)

        # Data
        obs_ = data.get_target_lightcurve(index).copy()
        if phase_window is not None:
            t0 = data.targets.data["t0"].loc[index]
            phase_window = np.asarray(phase_window)+t0
            obs_ = obs_[obs_["time"].between(*phase_window)]

        coef = 10 ** (-(obs_["zp"] - zp) / 2.5)
        obs_["flux_zp"] = obs_["flux"] * coef
        obs_["fluxerr_zp"] = obs_["fluxerr"] * coef

        # Model
        if bands is None:
            bands = np.unique(obs_["band"])

        # = axes and figure = #
        if ax is None:
            if fig is None:
                import matplotlib.pyplot as plt
                fig = plt.figure(figsize=[7,4])
            ax = fig.add_subplot(111)
        else:
            fig = ax.figure
        
        
        colors = get_band_color(bands)
        if show_truth:
            fig = my_show_lightcurve(data.targets, bands, ax=ax, fig=fig, index=index, 
                                           format_time=format_time, t0_format=t0_format, 
                                           zp=zp, colors=colors,
                                           zorder=2, 
                                           **lc_prop)
        elif format_time:
            from matplotlib import dates as mdates        
            locator = mdates.AutoDateLocator()
            formatter = mdates.ConciseDateFormatter(locator)
            ax.xaxis.set_major_locator(locator)
            ax.xaxis.set_major_formatter(formatter)
        else:
            ax.set_xlabel("time [in day]", fontsize="large")



        for band_, color_ in zip(bands, colors):
            if color_ is None:
                ecolor = to_rgba("0.4", 0.2)
            else:
                ecolor = to_rgba(color_, 0.2)
                
            obs_band = obs_[obs_["band"] == band_]
            times = obs_band["time"] if not format_time else Time(obs_band["time"], format=t0_format).datetime
            ax.scatter(times, obs_band["flux_zp"],
                       color=color_, zorder=4, s=3, **kwargs)
            ax.errorbar(times, obs_band["flux_zp"],
                        yerr= obs_band["fluxerr_zp"],
                        ls="None", marker="None", ecolor=ecolor, 
                        zorder=3,
                        **kwargs)
            # Set the x-limits to match the defined phase window
            if phase_window is not None:
                ax.set_xlim(Time(phase_window, format=t0_format).datetime)
            ax.set_ylim(bottom=10**(-3), top=None)
        return fig


def my_show_lightcurve(targets, band, index, params=None,
                            ax=None, fig=None, colors=None,
                            time_range=[-20,50], npoints=500,
                            zp=25, zpsys="ab",
                            format_time=True, t0_format="mjd", 
                            in_mag=False, invert_mag=True, **kwargs):
    if params is None:
        params = {}
            
    template = targets.get_target_template(index, **params)
    from skysurvey.config import get_band_color
    # get the sncosmo_model
    if params is None:
        params = {}

    sncosmo_model = template.get(**params)
        
    # ------- #
    #  x-data #
    # ------- #
    # time range
    t0 = sncosmo_model.get("t0")
    times = np.linspace(*np.asarray(time_range)+t0, npoints)

    # ------- #
    #  y-data #
    # ------- #        
    # flux
    band = np.atleast_1d(band)
    values = template.get_lightcurve(band,
                                 times, in_mag=in_mag,
                                 zp=zp, zpsys=zpsys,
                                 sncosmo_model=sncosmo_model)
        
    # ------- #
    #  axis   #
    # ------- #                    
    if ax is None:
        if fig is None:
            import matplotlib.pyplot as plt
            fig = plt.figure(figsize=[7,4])
        ax = fig.add_subplot(111)
    else:
        fig = ax.figure

    # ------- #
    #  Plot   #
    # ------- #  
    # The plot
    if format_time:
        from astropy.time import Time
        times = Time(times, format=t0_format).datetime

    colors = np.resize(colors, len(values))
    for band_, value_, color_ in zip(band, values, colors):
        if color_ is None: # default back to config color
            color_ = get_band_color(band_)

        ax.plot(times, value_, color=color_, **kwargs)

    # ------- #
    #  Format #
    # ------- #  
    # mag upside down
    if in_mag and invert_mag:
        ax.invert_yaxis()
    # time format
    if format_time:
        from matplotlib import dates as mdates        
        locator = mdates.AutoDateLocator()
        formatter = mdates.ConciseDateFormatter(locator)
        ax.xaxis.set_major_locator(locator)
        ax.xaxis.set_major_formatter(formatter)
    else:
        ax.set_xlabel("time [in day]", fontsize="large")

    if in_mag:
        ax.set_ylabel(f"Magnitude", fontsize="large")
    elif zp is None:
        ax.set_ylabel(f"Flux [erg/s/cm^2/A]", fontsize="large")
    else:
        ax.set_ylabel(f"Flux [zp={zp}]", fontsize="large") #
    ax.set_yscale('log')
    return fig


def my_mag_lightcurve(targets, band, index, params=None,
                            ax=None, fig=None, colors=None,
                            time_range=[-20,50], npoints=500,
                            zp=25, zpsys="ab",
                            format_time=True, t0_format="mjd", 
                            in_mag=False, invert_mag=True, **kwargs):
    if params is None:
        params = {}
            
    template = targets.get_target_template(index, **params)
    from skysurvey.config import get_band_color
    # get the sncosmo_model
    if params is None:
        params = {}

    sncosmo_model = template.get(**params)
        
    # ------- #
    #  x-data #
    # ------- #
    # time range
    t0 = sncosmo_model.get("t0")
    times = np.linspace(*np.asarray(time_range)+t0, npoints)

    # ------- #
    #  y-data #
    # ------- #        
    # flux
    band = np.atleast_1d(band)
    values = template.get_lightcurve(band,
                                 times, in_mag=in_mag,
                                 zp=zp, zpsys=zpsys,
                                 sncosmo_model=sncosmo_model)
        
    # ------- #
    #  axis   #
    # ------- #                    
    if ax is None:
        if fig is None:
            import matplotlib.pyplot as plt
            fig = plt.figure(figsize=[7,4])
        ax = fig.add_subplot(111)
    else:
        fig = ax.figure

    # ------- #
    #  Plot   #
    # ------- #  
    # The plot
    if format_time:
        from astropy.time import Time
        times = Time(times, format=t0_format).datetime

    colors = np.resize(colors, len(values))
    for band_, value_, color_ in zip(band, values, colors):
        if color_ is None: # default back to config color
            color_ = get_band_color(band_)

        ax.plot(times, value_, color=color_, **kwargs)

    # ------- #
    #  Format #
    # ------- #  
    # mag upside down
    if in_mag and invert_mag:
        ax.invert_yaxis()
    # time format
    if format_time:
        from matplotlib import dates as mdates        
        locator = mdates.AutoDateLocator()
        formatter = mdates.ConciseDateFormatter(locator)
        ax.xaxis.set_major_locator(locator)
        ax.xaxis.set_major_formatter(formatter)
    else:
        ax.set_xlabel("time [in day]", fontsize="large")

    if in_mag:
        ax.set_ylabel(f"Magnitude", fontsize="large")
    elif zp is None:
        ax.set_ylabel(f"Flux [erg/s/cm^2/A]", fontsize="large")
    else:
        ax.set_ylabel(f"Flux [zp={zp}]", fontsize="large")
    ax.plot(times, [22.4]*len(times), color='orange', label='ULTRASAT sensitivity')
    ax.plot(times, [24]*len(times), color='black', label='LSST sensitivity')
    ax.legend()
    return fig
