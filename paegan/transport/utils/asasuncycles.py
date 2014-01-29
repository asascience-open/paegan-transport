from datetime import timedelta
import numpy as np
import pytz


class SunCycles(object):

    SETTING = "sunset"
    RISING = "sunrise"

    @classmethod
    def cycles(cls, **kwargs):
        """
        Classmethod for convienence in returning both the sunrise and sunset
        based on a location and date.  Always calculates the sunrise and sunset on the
        given date, no matter the time passed into the function in the datetime object.

        Parameters:
            loc = Location4D (object)
                OR
            point = Shapely point (object)
            time = datetime in UTC (object)
                OR
            lat = latitude (float)
            lon = longitude (float)
            time = datetime in UTC (object)

        Returns:
            { 'sunrise': datetime in UTC, 'sunset': datetime in UTC }


        Sources:
            http://williams.best.vwh.net/sunrise_sunset_example.htm
        """
        if "loc" not in kwargs:
            if "point" not in kwargs:
                if "lat" not in kwargs or "lon" not in kwargs:
                    raise ValueError("You must supply some form of lat/lon coordinates")
                else:
                    lat = kwargs.get("lat")
                    lon = kwargs.get("lon")
            else:
                lat = kwargs.get("point").y
                lon = kwargs.get("point").x
            if "time" not in kwargs:
                raise ValueError("You must supply a datetime object")
            else:
                time = kwargs.get("time")
        else:
            lat = kwargs.get("loc").latitude
            lon = kwargs.get("loc").longitude
            time = kwargs.get("loc").time

        # Convert time to UTC.  Save passed in timezone to return later.
        if time.tzinfo is None:
            time = time.replace(tzinfo=pytz.utc)
            original_zone = pytz.utc
        else:
            original_zone = time.tzinfo

        local_jd = time.timetuple().tm_yday
        utc_jd = time.astimezone(pytz.utc).timetuple().tm_yday

        # We ALWAYS want to return the sunrise/sunset for the day that was passed
        # in (with timezone accounted for), regardless of what the UTC day is.  Modify
        # the UTC julian day here if need be.
        comp = cmp(utc_jd, local_jd)
        if comp == 1:
            utc_jd -= 1
        elif comp == -1:
            utc_jd += 1

        time = time.replace(hour=0, minute=0, second=0, microsecond=0)

        rising_h, rising_m = cls._calc(jd=utc_jd, lat=lat, lon=lon, stage=cls.RISING)
        setting_h, setting_m = cls._calc(jd=utc_jd, lat=lat, lon=lon, stage=cls.SETTING)

        # _calc returns UTC hours and minutes, so assume time is in UTC for a few lines...
        rising = time.replace(tzinfo=pytz.utc) + timedelta(hours=rising_h, minutes=rising_m)
        setting = time.replace(tzinfo=pytz.utc) + timedelta(hours=setting_h, minutes=setting_m)
        # LOOK: We may be adding 24 hours to the setting time.  Why?
        if setting < rising:
            setting = setting + timedelta(hours=24)

        rising = rising.astimezone(original_zone)
        setting = setting.astimezone(original_zone)

        return { cls.RISING : rising, cls.SETTING : setting}

    @classmethod
    def _calc(cls, **kwargs):
        """
        Calculate sunrise or sunset based on:
        Parameters:
            jd: Julian Day
            lat: latitude
            lon: longitude
            stage: sunrise or sunset
        """
        zenith = 90.833333 # offical value

        jd = kwargs.get("jd", None)
        lat = kwargs.get("lat", None)
        lon = kwargs.get("lon", None)
        stage = kwargs.get("stage", None)
        if jd is None or stage is None or lat is None or lon is None:
            raise ValueError("Must supply an 'jd', 'lat, 'lon', and 'stage' parameter")

        if stage != SunCycles.RISING and stage != SunCycles.SETTING:
            raise ValueError("'stage' parameter must be %s or %s" % (SunCycles.RISING, SunCycles.SETTING))

        longhr = lon / 15.
        if stage == SunCycles.RISING:
            apx = jd + ( (6 - longhr) / 24 )
        elif stage == SunCycles.SETTING:
            apx = jd + ( (18 - longhr) / 24 )

        sun_mean_anom = ( 0.9856 * apx ) - 3.289 # sun's mean anomaly
        #sun's longitude
        sun_lon = sun_mean_anom + (1.916 * np.sin( np.radians(sun_mean_anom) )) \
                + (0.02 * np.sin( np.radians(2 * sun_mean_anom) )) + 282.634

        if sun_lon > 360:
            sun_lon = sun_lon - 360
        elif sun_lon < 0:
            sun_lon = sun_lon + 360

        right_ascension = np.degrees(np.arctan( 0.91764 * np.tan( np.radians(sun_lon) ) )) # sun's right ascension

        if right_ascension > 360:
            right_ascension = right_ascension - 360
        elif right_ascension < 0:
            right_ascension = right_ascension + 360

        # put sun's right ascension value in the same quadrant as the sun's
        # true longitude
        lQuad = 90. * np.floor(sun_lon / 90.)
        raQuad = 90. * np.floor(right_ascension / 90.)
        right_ascension = right_ascension + ( lQuad - raQuad)
        right_ascension = right_ascension / 15. # Convert to hours

        # Sun's declination
        sinDecl = 0.39782 * np.sin( np.radians(sun_lon) )
        cosDecl = np.cos( np.arcsin( sinDecl ) )

        # Sun's local hour angle
        cosHr = (np.cos( np.radians(zenith) ) - ( sinDecl * np.sin(np.radians(lat)) )) \
                / ( cosDecl * np.cos( np.radians(lat) ) )

        if cosHr > 1: # Sun doesnt rise on this loc on this date
            return -1, -1
        elif cosHr < -1: # Sun doesnt set on this location on this date
            return -1, -1
        elif stage == SunCycles.RISING: # Sunrise
            hr = 360 - np.degrees(np.arccos(cosHr))
        elif stage == SunCycles.SETTING:  # Sunset
            hr = np.degrees(np.arccos(cosHr))

        hr = hr / 15. # Convert angle to hours

        localTime = hr + right_ascension - ( 0.06571 * apx ) - 6.622# local meantime of rise/set

        UTtime = localTime - longhr # adjust to UTC

        if UTtime < 0:
            UTtime = UTtime + 24
        elif UTtime > 24:
            UTtime = UTtime - 24

        hour = np.floor(UTtime)

        minute = (UTtime - hour) * 60
        if minute == 60:
            hour = hour + 1
            minute = 0

        return hour, minute
