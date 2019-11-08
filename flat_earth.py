def flat_earth(lon0, lat0, bearing=None, distance=None,
               distance_east=None, distance_north=None,
               lat1=None, lon1=None):

    """Flat earth approximation: From reference position (lon0, lat0) and
    a second positon in the vicinity thereof to distance and bearing
    (and vice versa).

    Approximation fails in the vicinity of either
    pole and at large distances.

    Fractional errors are of order (distance/R)**2.

    To compare this solution with an excact one from the geographiclib
    module, execute the following lines of Python code:
    import math
    from geographiclib.geodesic import Geodesic
    geod = Geodesic.WGS84
    g = geod.Inverse(lat0, lon0, lat1, lon1)
    print "The distance is {:.3f} m.".format(g['s12'])
    print "The azimuth is {:.3f} m.".format((450.-g['azi1']) % 360.)
    g = geod.Direct(lat0, lon0, bearing, distance)
    print "The position is ({:.8f}, {:.8f}).".format(g['lat2'],g['lon2'])
    
    Source for the following equations:
    http://edwilliams.org/avform.htm (linked to by
    https://www.nhc.noaa.gov)

    Args:

    lon0 (float) - Reference longitude in decimal degrees

    lat0 (float) - Reference latitude in decimal degrees

    bearing (float, optional) - Bearing with reference position as
                                origin in deg clockwise from north;
                                default is None

    distance (float, optional) - Distance from reference position in
                                 meters; default is None

    distance_east(float, optional) - Distance along east-west axis
                                     from reference position in
                                     meters; default is None

    distance_north(float, optional) - Distance along north-south axis
                                      from reference position in
                                      meters; default is None

    lon1 (float, optional) - Longitude of second position in vicinity
                             to the reference position in decimal
                             degrees; default is None

    lat1 (float, optional) - Latitude of second position in vicinity
                             to the reference position in decimal
                             degrees; default is None

    """

    import numpy as np

    f = 1./298.257223563 # Flattening (WGS 84)
    a = 6378137. # [m] Equatorial radius of Earth (WGS 84)
    e_squared = f*(2-f)
    lat0_rad = np.radians(lat0)
    r1 = (a*(1.-e_squared)) / (1.-e_squared*np.sin(lat0_rad)**2)**(3./2.) # Meridional radius of curvature
    r2 = a / np.sqrt(1.-e_squared*np.sin(lat0_rad)**2) # Radius of curvarture in the prime vertical

    class ReturnValue(object):
        pass
        
    def inverse(lon1, lat1):
        d_lat = lat1-lat0
        d_lon = lon1-lon0
        dist_n = r1*np.radians(d_lat)
        dist_e = r2*np.cos(lat0_rad)*np.radians(d_lon)
        dist = np.sqrt(dist_n**2 + dist_e**2)
        bearing = np.degrees(np.arctan2(dist_n, dist_e))
        result = ReturnValue()
        result.dist_n = dist_n
        result.dist_e = dist_e
        result.dist = dist
        result.bearing = bearing
        return result
        
    def direct(dist_e, dist_n):
        d_lat_rad = dist_n/r1
        d_lon_rad = dist_e/(r2*np.cos(lat0_rad))
        result = ReturnValue()
        result.lat = lat0+np.degrees(d_lat_rad)
        result.lon = lon0+np.degrees(d_lon_rad)
        return result

    if lon1 is not None and lat1 is not None:
        return inverse(lon1, lat1)

    if distance_east is not None and distance_north is not None:
        return direct(distance_east, distance_north)
       
    if bearing is not None and distance is not None:
        distance_east = distance*np.cos(np.radians(bearing))
        distance_north = distance*np.sin(np.radians(bearing))
        return direct(distance_east, distance_north)
