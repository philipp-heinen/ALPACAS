import numpy as np
import BasicFunctions as bf
import Plotting as pl
import CelestialMechanics as cm
import Constants as const

class Sundial:
    def __init__(self, latitude=45, longitude=0, orientation="horizontal"):
        if -90 <= latitude <= 90:
            self.latitude = latitude
        else:
            raise Exception("latitude must be between -90 and 90.")
            
        if -180 <= longitude <= 180:
            self.latitude = latitude
        else:
            raise Exception("longitude must be between -180 and 180.")
     
        if isinstance(orientation, str) and orientation in ["horizontal", "vertical", "equatorial_upper_side", "equatorial_lower_side", "polar"]:
            self.orientation = orientation
        elif isinstance(orientation, list) and len(orientation) == 2:
            if -90 <= orientation[0] <= 90 and -90 <= orientation[1] <= 90:
                self.orientation = orientation
            else: 
                raise Exception("Both orientation angles must be between -90 and 90")
        else:
            raise Exception("orientation must be either one out of 'horizontal', 'vertical', 'equatorial_upper_side', 'equatorial_lower_side', 'polar' or a list of length 2 with the angles of the dial.")
         
        self.phi = latitude/180*np.pi
        self.chi = longitude/180*np.pi
        
        if self.phi >= 0:
            self.nh = True
        else:
            self.nh = False
         
        #this value determines how precisely the end point(s) of straight lines that are not defined for all declinations are determined
        self.eps = 0.01

        
    def _nodus(self, tau, delta):
        if self.orientation == "horizontal":
            return bf.nodus_horizontal(tau, delta, self.phi)
        elif self.orientation == "vertical":
            return bf.nodus_vertical(tau, delta, self.phi)
        elif self.orientation == "equatorial_upper_side":
            return bf.nodus_equatorial(tau, delta, self.phi, upper_side=True)
        elif self.orientation == "equatorial_lower_side":
            return bf.nodus_equatorial(tau, delta, self.phi, upper_side=False)
        elif self.orientation == "polar":
            return bf.nodus_polar(tau, delta, self.phi)
        else:
            return bf.nodus_arbitrary_orientation(tau, delta, self.phi, self.orientation[0], self.orientation[1])

    def _gnomon(self, tau):
        if self.orientation == "horizontal":
            return bf.gnomon_horizontal(tau, self.phi)
        elif self.orientation == "vertical":
            return bf.gnomon_vertical(tau, self.phi)
        elif self.orientation == "equatorial_upper_side":
            return bf.gnomon_equatorial(tau, self.phi, upper_side=True)
        elif self.orientation == "equatorial_lower_side":
            return bf.gnomon_equatorial(tau, self.phi, upper_side=False)
        elif self.orientation == "polar":
            raise Exception("The polar gnomon is parallel to the dial of a polar sundial.")
        else:
            return bf.gnomon_arbitrary_orientation(tau, self.phi, self.orientation[0], self.orientation[1])
    
    def _nodus_vec(self, taus, delta):
        X = np.array([])
        Y = np.array([])
        for tau in taus:
            nd = self._nodus(tau, delta)
            X = np.append(X, nd[0])
            Y = np.append(Y, nd[1])
        return X, Y
    
    def _declination_vec(self, ts):
        delta = np.array([])
        for t in ts:
            delta = np.append(delta, cm.declination(t))
        return delta
    
    def _eot_vec(self, ts):
        eot = np.array([])
        for t in ts:
            eot = np.append(eot, cm.declination(t))
        return eot
    
    def _origin_polar_gnomon(self):
        if self.orientation == "horizontal":
            opg = bf.position_gnomon_horizontal(self.phi)
        elif self.orientation == "vertical":
            opg = bf.position_gnomon_vertical(self.phi)
        elif self.orientation == "equatorial_upper_side" or self.orientation == "equatorial_lower_side":
            opg = np.array([0,0])
        elif self.orientation == "polar":
            raise Exception("The polar gnomon is parallel to the dial of a polar sundial.")
        else:
            opg = bf.position_gnomon_arbitrary_orientation(self.phi, self.orientation[0], self.orientation[1])
        return opg

    def _max_decs(self, f):
        f1, f2 = (np.array([]), np.array([]))
        decs = np.arange(-const.ECLIPTIC/180*np.pi, const.ECLIPTIC/180*np.pi+self.eps, self.eps)
        if not self.nh:
            decs = np.flip(decs)
        for dec in decs:
            tmp = f(dec)
            f1 = np.append(f1, tmp[0])
            f2 = np.append(f2, tmp[1])
        finite_ind = np.where(np.isfinite(f1)&np.isfinite(f2))[0]
        if len(finite_ind) < 2:
            return np.nan, np.nan
        else:
            return decs[finite_ind[0]], decs[finite_ind[-1]]
          
    def init_dial_plot(self, xsize=10, ysize=10, cm_per_unit=2):
        """
        Initiates a new dial plot.

        Parameters
        ----------
        xsize : the extension of the plot in x direction, in units of the (nodal) gnomon length.
        ysize : the extension of the plot in y direction, in units of the (nodal) gnomon length.
        cm_per_unit : centimeters per unit length.

        Returns
        -------
        None.

        """
        self.plotter = pl.Plotter(xsize, ysize, cm_per_unit)
        self.xsize = xsize
        self.ysize = ysize
        
    def save_dial_plot(self, filename, resolution=300):
        """
        Saves a dial plot to a file 'filename' with a resolution provided in dpi. 
        """
        self.plotter.save(filename=filename, resolution=resolution)
        
    def add_text(self, text, x, y, fontsize=12, color="black"):
        """
        Adds text at position (x,y). 
        """
        self.plotter.text(text, x, y, fontsize=fontsize, color=color)
        
    def add_nodus_pos(self, style=None):
        """
        Adds a marker for the position of the nodal gnomon.
        """
        self.plotter.plot_point([0,0], style=style)
        
    def add_origin_polar_gnomon(self, style=None):
        """
        Adds a marker for the position of the polar gnomon.
        """
        self.plotter.plot_point(self._origin_polar_gnomon(), style)
        
    def add_date_line(self, date, style=None, label=None, labelstyle=None, rendering=0.01):
        """
        Draws a date line for the specified date.

        Parameters
        ----------
        date : Can be either "equinox", "summer_solstice", "winter_solstice" or a custom date provided as [month, day].
        style : line style, provided as dictionary.
        label : the line label.
        labelstyle: style of the label, provided as dictionary.
        rendering : the rendering of the line, expressed as hour angle.

        Returns
        -------
        None.

        """
        if date == "equinox":
            delta = 0
        elif date == "summer_solstice":
            if self.nh:
                delta = const.ECLIPTIC/180*np.pi
            else:
                delta = -const.ECLIPTIC/180*np.pi
        elif date == "winter_solstice":
            if self.nh:
                delta = -const.ECLIPTIC/180*np.pi
            else:
                delta = const.ECLIPTIC/180*np.pi
        elif isinstance(date, (list, np.ndarray)) and len(date) == 2:
            delta = cm.declination_date(date[0], date[1])
        else:
            raise Exception("Date must be either 'equinox', 'summer_solstice', 'winter_solstice' or a date specified as [month, day]")
            
        taumax = bf.sunset(delta, self.phi)
     
        taus = np.append(np.flip(-np.arange(rendering, taumax, rendering)), np.arange(0, taumax, rendering))
        
        if labelstyle is None:
            labelstyle = {}
        if self.nh: 
            labelstyle.setdefault("rotation_type", "parallel_up")
            labelstyle.setdefault("side", "left")
        else:
            labelstyle.setdefault("rotation_type", "parallel_down")
            labelstyle.setdefault("side", "right")
        labelstyle.setdefault("tightness", 0.15)
        labelstyle.setdefault("position", 0.8)
        
        X, Y = self._nodus_vec(taus, delta)
        self.plotter.plot_curve(X, Y, style=style, label=label, labelstyle=labelstyle)
           
    def _get_hours_and_labels(self, which, which_labels, start_at_zero=False):
        if which == "hourly":
            hours = np.arange(0, 24, 1)
            if which_labels == "roman_numerals":
                labels = [pl.roman_num[hour] for hour in hours]
            elif which_labels == "arabic_numerals":
                labels = [str(hour) for hour in hours]
            elif which_labels == "no_labels":
                labels = [None for hour in hours]
            else:
                raise Exception("In combination with which='hourly', the labeltype must be 'roman_numerals', 'arabic_numerals' or 'no_labels'.")
        elif which == "half-hourly":
            hours = np.arange(0, 24, 0.5)
            labels = []
            if which_labels == "roman_numerals":
                for i in range(24):
                    labels.append(pl.roman_num[i])
                    labels.append(None)
            elif which_labels == "arabic_numerals":
                for i in range(24):
                    labels.append(str(i))
                    labels.append(None)
            elif which_labels == "no_labels":
                labels = [None for i in range(48)]
            else:
                raise Exception("In combination with which='half-hourly', the labeltype must be 'roman_numerals', 'arabic_numerals' or 'no_labels'.")                
        elif which == "two-hourly":
            hours = np.arange(0, 24, 2)
            if which_labels == "roman_numerals":
                labels = [pl.roman_num[hour] for hour in hours]
            elif which_labels == "arabic_numerals":
                labels = [str(hour) for hour in hours]
            elif which_labels == "no_labels":
                labels = [None for hour in hours]
            else:
                raise Exception("In combination with which='two-hourly', the labeltype must be 'roman_numerals', 'arabic_numerals' or 'no_labels'.")
        elif isinstance(which, (list, np.ndarray)):
            hours = np.asarray(which)
            if which_labels == "no_labels":
                labels = [None for hour in hours]
            elif isinstance(which_labels, (list, np.ndarray)) and len(which_labels) == len(which):
                labels = which_labels
            else:
                raise Exception("If a list is provided for which, which_labels must be either a list of the same length containing the labels or 'no_labels'. ")
        else:
            raise Exception("which must be 'hourly', 'half-hourly', 'two-hourly' or a list with the desired hours")        
        
        if not start_at_zero:
            hours -= 12
        
        return hours, np.asarray(labels)
        
    def _add_apparent_solar_time(self, taus, labels, gnomon_type="nodus", style=None, labelstyle=None):
        if self.nh:
            decl_summer_solstice = const.ECLIPTIC/180*np.pi
            decl_winter_solstice = -const.ECLIPTIC/180*np.pi
        else:
            decl_summer_solstice = -const.ECLIPTIC/180*np.pi
            decl_winter_solstice = const.ECLIPTIC/180*np.pi

        maxtau_max = bf.sunset(decl_summer_solstice, self.phi)
        maxtau_min = bf.sunset(decl_winter_solstice, self.phi)

        labels = np.delete(labels, np.where((taus>maxtau_max)|(taus<-maxtau_max)))
        taus = np.delete(taus, np.where((taus>maxtau_max)|(taus<-maxtau_max)))
        
        if gnomon_type == "nodus":
            for tau, label in zip(taus, labels):
                if np.all(np.isfinite(self._nodus(tau, decl_summer_solstice))) and np.all(np.isfinite(self._nodus(tau, decl_winter_solstice))):
                    self.plotter.plot_straight_line(self._nodus(tau, decl_summer_solstice), self._nodus(tau, decl_winter_solstice), style=style, label=label, labelstyle=labelstyle)
                else:
                    decmin, decmax = self._max_decs(lambda x: self._nodus(tau, x))
                    self.plotter.plot_straight_line(self._nodus(tau, decmax), self._nodus(tau, decmin), style=style, label=label, labelstyle=labelstyle)
        elif gnomon_type == "polar":
            if self.orientation == "polar":
                raise Exception("For a polar sundial, gnomon_type must be 'nodus'.")
            for tau, label in zip(taus, labels):
                    length = max(self.xsize, self.ysize)*2
                    angle = self._gnomon(tau)
                    origin = self._origin_polar_gnomon()
                    self.plotter.plot_straight_line(origin, origin+length*np.array([np.sin(angle), np.cos(angle)]), style=style, label=label, labelstyle=labelstyle)
        else:
            raise Exception("gnomon_type must be 'nodus' or 'polar'")
                                             
    def add_apparent_solar_time(self, which="hourly", gnomon_type="nodus", style=None, which_labels="roman_numerals", labelstyle=None):
        """
        Adds lines indicating the apparent solar time.

        Parameters
        ----------
        which : Which hours to draw. Can be "hourly", "half-hourly", "two-hourly" or a costum list of hours. Default is "hourly".
        gnomon_type : Indicates whether to draw lines for a "polar" or "nodus" gnomon. Default is "nodus".
        style : line style, provided as dictionary.
        which_labels : Which labels to assign to the hour lines. Can be "no_labels", "roman_numerals", "arabic_numerals" or a custom list of strings whose length must match that of which.
        labelstyle : label style, provided as dictionary.

        Returns
        -------
        None.

        """
        hours, labels = self._get_hours_and_labels(which, which_labels)
        taus = hours*np.pi/12
        self._add_apparent_solar_time(taus=taus, labels=labels, gnomon_type=gnomon_type, style=style, labelstyle=labelstyle)           

    def add_apparent_zone_time(self, which="hourly", timezone="closest", gnomon_type="nodus", style=None, which_labels="roman_numerals", labelstyle=None):
        """
        Adds lines indicating the apparent zone time.

        Parameters
        ----------
        which : Which hours to draw. Can be "hourly", "half-hourly", "two-hourly" or a costum list of hours. Default is "hourly".
        timezone: Can be either "closest", i.e. the time zone is determined based on the longitude, or as an integer.
        gnomon_type : Indicates whether to draw lines for a "polar" or "nodus" gnomon. Default is "nodus".
        style : line style, provided as dictionary.
        which_labels : Which labels to assign to the hour lines. Can be "no_labels", "roman_numerals", "arabic_numerals" or a custom list of strings whose length must match that of which.
        labelstyle : label style, provided as dictionary.

        Returns
        -------
        None.

        """        
        hours, labels = self._get_hours_and_labels(which, which_labels)
        taus = hours*np.pi/12
        if timezone == "closest":
            offset = self.chi-np.round(self.chi/(np.pi/12))*(np.pi/12)
        else:
            offset = self.chi-timezone*(np.pi/12)
        taus += offset
        self._add_apparent_solar_time(taus=taus, labels=labels, gnomon_type=gnomon_type, style=style, labelstyle=labelstyle)   

    def _add_mean_solar_time(self, taus, labels, half_year="spring", style=None, labelstyle=None, rendering=0.01):
        if self.nh:
            decl_summer_solstice = const.ECLIPTIC/180*np.pi
            decl_winter_solstice = -const.ECLIPTIC/180*np.pi
        else:
            decl_summer_solstice = -const.ECLIPTIC/180*np.pi
            decl_winter_solstice = const.ECLIPTIC/180*np.pi

        maxtau_max = bf.sunset(decl_summer_solstice, self.phi)
        maxtau_min = bf.sunset(decl_winter_solstice, self.phi)

        labels = np.delete(labels, np.where((taus>maxtau_max)|(taus<-maxtau_max)))
        taus = np.delete(taus, np.where((taus>maxtau_max)|(taus<-maxtau_max)))   
        
        if half_year == "spring":
            ts = np.flip(np.append(np.arange(cm.find_solstice("winter"), 1, rendering), np.arange(0, cm.find_solstice("summer"), rendering)))
        elif half_year == "autumn":
            ts = np.arange(cm.find_solstice("summer"), cm.find_solstice("winter"), rendering)

        for tau, label in zip(taus, labels):
            X, Y = (np.array([]), np.array([]))
            for t in ts: 
                tmp = self._nodus(tau+cm.equation_of_time(t)*np.pi/720, cm.declination(t))
                X = np.append(X, tmp[0])
                Y = np.append(Y, tmp[1])
            self.plotter.plot_curve(X, Y, style=style, label=label, labelstyle=labelstyle)
        
    def add_mean_solar_time(self, which="hourly", half_year="spring", style=None, which_labels="roman_numerals", labelstyle=None, rendering=0.01):
        """
        Adds lines indicating the mean solar time.

        Parameters
        ----------
        which : Which hours to draw. Can be "hourly", "half-hourly", "two-hourly" or a costum list of hours. Default is "hourly".
        half_year: specifies whether lines should be valid for "spring" half-year (between winter and summer solstice) or "autumn" half-year (between summer and winter solstice).
        style : line style, provided as dictionary.
        which_labels : Which labels to assign to the hour lines. Can be "no_labels", "roman_numerals", "arabic_numerals" or a custom list of strings whose length must match that of which.
        labelstyle : label style, provided as dictionary.
        rendering: the rendering of the line, expressed as hour angle.

        Returns
        -------
        None.

        """
        hours, labels = self._get_hours_and_labels(which, which_labels)
        taus = hours*np.pi/12
        self._add_mean_solar_time(taus=taus, labels=labels, half_year=half_year, style=style, labelstyle=labelstyle, rendering=rendering) 
        
    def add_mean_zonal_time(self, which="hourly", timezone="closest", half_year="spring", style=None, which_labels="roman_numerals", labelstyle=None, rendering=0.01):
        """
        Adds lines indicating the mean zonal (=civil) time.

        Parameters
        ----------
        which : Which hours to draw. Can be "hourly", "half-hourly", "two-hourly" or a costum list of hours. Default is "hourly".
        timezone: Can be either "closest", i.e. the time zone is determined based on the longitude, or as an integer.
        half_year: specifies whether lines should be valid for "spring" half-year (between winter and summer solstice) or "autumn" half-year (between summer and winter solstice).
        style : line style, provided as dictionary.
        which_labels : Which labels to assign to the hour lines. Can be "no_labels", "roman_numerals", "arabic_numerals" or a custom list of strings whose length must match that of which.
        labelstyle : label style, provided as dictionary.
        rendering: the rendering of the line, expressed as hour angle.

        Returns
        -------
        None.

        """        
        hours, labels = self._get_hours_and_labels(which, which_labels)
        taus = hours*np.pi/12
        if timezone == "closest":
            offset = self.chi-np.round(self.chi/(np.pi/12))*(np.pi/12)
        else:
            offset = self.chi-timezone*(np.pi/12)
        taus += offset
        self._add_mean_solar_time(taus=taus, labels=labels, half_year=half_year, style=style, labelstyle=labelstyle, rendering=rendering)
        
    def add_babylonian_hours(self, which="hourly", style=None, which_labels="roman_numerals", labelstyle=None):
        """
        Adds lines indicating Babylonian hours.

        Parameters
        ----------
        which : Which hours to draw. Can be "hourly", "half-hourly", "two-hourly" or a costum list of hours. Default is "hourly".
        style : line style, provided as dictionary.
        which_labels : Which labels to assign to the hour lines. Can be "no_labels", "roman_numerals", "arabic_numerals" or a custom list of strings whose length must match that of which.
        labelstyle : label style, provided as dictionary.

        Returns
        -------
        None.

        """
        if (self.phi > np.pi/2-const.ECLIPTIC/180*np.pi) or (self.phi < -np.pi/2+const.ECLIPTIC/180*np.pi):
            raise Exception("Babylonian hours are not supported for latitudes beyond the polar circles.")
        
        hours, labels = self._get_hours_and_labels(which, which_labels, start_at_zero=True)
        
        if self.nh:
            decl_summer_solstice = const.ECLIPTIC/180*np.pi
            decl_winter_solstice = -const.ECLIPTIC/180*np.pi
        else:
            decl_summer_solstice = -const.ECLIPTIC/180*np.pi
            decl_winter_solstice = const.ECLIPTIC/180*np.pi

        maxtau_max = bf.sunset(decl_summer_solstice, self.phi)
        maxtau_min = bf.sunset(decl_winter_solstice, self.phi)
            
        for hour, label in zip(hours, labels):
            tau_summer = -maxtau_max+hour*np.pi/12
            tau_winter = -maxtau_min+hour*np.pi/12
            if -maxtau_max < tau_summer < maxtau_max:
                if np.all(np.isfinite(self._nodus(tau_summer, decl_summer_solstice))) and np.all(np.isfinite(self._nodus(tau_winter, decl_winter_solstice))):
                    self.plotter.plot_straight_line(self._nodus(tau_summer, decl_summer_solstice), self._nodus(tau_winter, decl_winter_solstice), style=style, label=label, labelstyle=labelstyle)
                else:    
                    decmin, decmax = self._max_decs(lambda x: self._nodus(-bf.sunset(x, self.phi)+hour*np.pi/12, x))
                    taumin, taumax = -bf.sunset(decmin, self.phi)+hour*np.pi/12, -bf.sunset(decmax, self.phi)+hour*np.pi/12
                    self.plotter.plot_straight_line(self._nodus(taumax, decmax), self._nodus(taumin, decmin), style=style, label=label, labelstyle=labelstyle)
                              
    def add_italian_hours(self, which="hourly", style=None, which_labels="roman_numerals", labelstyle=None):
        """
        Adds lines indicating Italian hours.

        Parameters
        ----------
        which : Which hours to draw. Can be "hourly", "half-hourly", "two-hourly" or a costum list of hours. Default is "hourly".
        style : line style, provided as dictionary.
        which_labels : Which labels to assign to the hour lines. Can be "no_labels", "roman_numerals", "arabic_numerals" or a custom list of strings whose length must match that of which.
        labelstyle : label style, provided as dictionary.

        Returns
        -------
        None.

        """
        if (self.phi > np.pi/2-const.ECLIPTIC/180*np.pi) or (self.phi < -np.pi/2+const.ECLIPTIC/180*np.pi):
            raise Exception("Italian hours are not supported for latitudes beyond the polar circles.")
        
        hours, labels = self._get_hours_and_labels(which, which_labels, start_at_zero=True)
        
        if self.nh:
            decl_summer_solstice = const.ECLIPTIC/180*np.pi
            decl_winter_solstice = -const.ECLIPTIC/180*np.pi
        else:
            decl_summer_solstice = -const.ECLIPTIC/180*np.pi
            decl_winter_solstice = const.ECLIPTIC/180*np.pi

        maxtau_max = bf.sunset(decl_summer_solstice, self.phi)
        maxtau_min = bf.sunset(decl_winter_solstice, self.phi)
            
        for hour, label in zip(hours, labels):
            tau_summer = maxtau_max-2*np.pi+hour*np.pi/12
            tau_winter = maxtau_min-2*np.pi+hour*np.pi/12
            if -maxtau_max < tau_summer < maxtau_max:
                if np.all(np.isfinite(self._nodus(tau_summer, decl_summer_solstice))) and np.all(np.isfinite(self._nodus(tau_winter, decl_winter_solstice))):
                    self.plotter.plot_straight_line(self._nodus(tau_summer, decl_summer_solstice), self._nodus(tau_winter, decl_winter_solstice), style=style, label=label, labelstyle=labelstyle)
                else:    
                    decmin, decmax = self._max_decs(lambda x: self._nodus(bf.sunset(x, self.phi)-2*np.pi+hour*np.pi/12, x))
                    taumin, taumax = bf.sunset(decmin, self.phi)-2*np.pi+hour*np.pi/12, bf.sunset(decmax, self.phi)-2*np.pi+hour*np.pi/12
                    self.plotter.plot_straight_line(self._nodus(taumax, decmax), self._nodus(taumin, decmin), style=style, label=label, labelstyle=labelstyle)
    
    def add_temporal_hours(self, which="hourly", style=None, which_labels="roman_numerals", labelstyle=None):
        """
        Adds lines indicating temporal hours.

        Parameters
        ----------
        which : Which hours to draw. Can be "hourly", "half-hourly", "two-hourly" or a costum list of hours. Default is "hourly".
        style : line style, provided as dictionary.
        which_labels : Which labels to assign to the hour lines. Can be "no_labels", "roman_numerals", "arabic_numerals" or a custom list of strings whose length must match that of which.
        labelstyle : label style, provided as dictionary.

        Returns
        -------
        None.

        """
        if (self.phi > np.pi/2-const.ECLIPTIC/180*np.pi) or (self.phi < -np.pi/2+const.ECLIPTIC/180*np.pi):
            raise Exception("Temporal hours are not supported for latitudes beyond the polar circles.")        
    
        hours, labels = self._get_hours_and_labels(which, which_labels, start_at_zero=True)

        labels = labels[np.where((hours>0)&(hours<12))]
        hours = hours[np.where((hours>0)&(hours<12))]
        
        if self.nh:
            decl_summer_solstice = const.ECLIPTIC/180*np.pi
            decl_winter_solstice = -const.ECLIPTIC/180*np.pi
        else:
            decl_summer_solstice = -const.ECLIPTIC/180*np.pi
            decl_winter_solstice = const.ECLIPTIC/180*np.pi

        maxtau_max = bf.sunset(decl_summer_solstice, self.phi)
        maxtau_min = bf.sunset(decl_winter_solstice, self.phi)
            
        for hour, label in zip(hours, labels):
            tau_summer = maxtau_max*(hour/6-1)
            tau_winter = maxtau_min*(hour/6-1)
            if -maxtau_max < tau_summer < maxtau_max:
                if np.all(np.isfinite(self._nodus(tau_summer, decl_summer_solstice))) and np.all(np.isfinite(self._nodus(tau_winter, decl_winter_solstice))):
                    self.plotter.plot_straight_line(self._nodus(tau_summer, decl_summer_solstice), self._nodus(tau_winter, decl_winter_solstice), style=style, label=label, labelstyle=labelstyle)
                else:    
                    decmin, decmax = self._max_decs(lambda x: self._nodus(bf.sunset(x, self.phi)*(hour/6-1), x))
                    taumin, taumax = bf.sunset(decmin, self.phi)*(hour/6-1), bf.sunset(decmax, self.phi)*(hour/6-1)
                    self.plotter.plot_straight_line(self._nodus(taumax, decmax), self._nodus(taumin, decmin), style=style, label=label, labelstyle=labelstyle)    
