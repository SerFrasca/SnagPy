import astropy
import time
import GD

### TIME ### ---------------------------

# mjd time ----------------------------

def mjd_now():
    pass


def mjd2v(mjd):
    pass


def v2mjd(v):
    pass


def diff_mjd(mjd1,mjd2):
    pass


def tt2mjd(tai):
    pass


def mjd2tt(mjd):
    pass


def mjd2s(mjd,typ):
    pass


def s2mjd(str):
    pass


def mjdyear(year):
    pass


def tai2mjd(tai):
    pass


def mjd2tai(mjd):
    pass



def adds2mjd(mjd1,nsec):
    pass



# gps time ----------------------------


def gps_now():
# Return the current GPS time as a float using Astropy.

    from astropy.time import Time

    return float(Time.now().gps)


def gps2mjd(tgps): # conversion from gps time to mjd
    pass


def gps2utc(gps):
    pass


def mjd2gps(mjd):
    pass


def v2gps():
    pass


# sidereal time

def sidereal_hour(t):
    pass


def gmst(t):
# Greenwich mean sidereal time (in hours)
    st=sidereal_hour(t)

    return st


def loc_sidtim(mjd,long):
    pass


def nearest_tsid(ts,t0):
    pass


def sidsolasid(t):
    pass


def nearest_tsid(ts,t0):
    pass



# other -----------------------------

def leap_seconds(mjd):
    pass


def tdt2tdb(mjd):
    pass


def tdt2tdb_hp(mjd):
    pass


def vday2doy(vday):
    pass


### ASTRO ###  __________________________________________________

# Coordinates --------------------------

def astro_coord(cin,cout,ai,di,ltsid,lat):
    pass


def astro_angle(in1,in2):
    pass


def spherical_distance(a1,d1,a2,d2):
    pass


def astro2rect(a,icrad):
    pass


def hour2deg(hour):
    pass



# Doppler --------------------------------

def doppler(time,fr0,alpha,delta,long,lat,coord):
    pass


def doppler2(time,fr0,alpha,delta,long,lat,coord):
    pass


def earth_v0(time,long,lat):
    pass


def earth_v1(time,long,lat):
    pass


def earth_v1_norot(time,coord):
    pass


def earth_v2(time,long,lat):
    pass


def v_detector(doptab,t):
    pass


def gw_doppler(doptab,source,t):
    pass


def gw_doppler1(v,source):
    pass


def dopp_from_gd(gdin,source):
    pass


def pos_detector(doptab,t):
    pass


def source_delay(doptab,sour,t):
    pass


def read_doppler_table(file,subsamp,fileout):
    pass


def reduce_doptab(doptab,tmin,tmax):
    pass



# Relativistic corrections ------------------------

def einst_effect(mjd,icplot):
    pass



def shapiro_delay(mjd,source):
    pass


# Kepler equation --------------------------

def kepler_elliptic_ofek(T,Q,Ecc):
    pass


def base_elliptic_orbit(ecc,sma,N):
    pass


# Other --------------------------------------

def tide(t1,lat,long,height):
    pass



def t_culm(long,rasc,t0):
    pass


def pss_astra(typ):
    pass


def galac_disc(n):
    pass


