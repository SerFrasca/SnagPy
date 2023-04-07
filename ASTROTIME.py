    # Copyright (C) 2023  Riccardo Felicetti, Sergio Frasca
    #  under GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007

'''
        Module ASTROTIME
Time and astronomy management

Sections:
>> TIME
> Basic time
> AstroPy time
> mjd time
> gps time
> sidereal time
>> ASTRO
> Constants
> Coordinates
> AstroPy coordinates
> Using Skyfield
> Using jpl tables
> Doppler
> Relativistic corrections
> Kepler equation
> Other
'''

import numpy as np
import astropy
from astropy import constants as const
from astropy.time import Time
from astropy import units as u
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.coordinates import Angle
from skyfield.api import load
import time
import datetime
import copy
import GD,BASIC,GWDATA

Pi=np.pi
deg2rad=Pi/180
rad2deg=1/deg2rad

### TIME ### ---------------------------

# basic time ---------------------------

def now(vt=0,typ='UTC',form='string'):
    '''
    Time now
    vt      vect time, if absent, present time
            vt=(year, month [1-12], day [1-31], hour, min, s, 0,0,0)
    typ    'UTC' or 'local'
    form   'string', 'struct', 'mjd', 'gps', 'unix'
    '''

    tz=time.timezone
    if vt == 0:
        t=time.time()
        vt=time.gmtime(t)
        vt=timestruct2v(vt)
    else:
        l=len(vt)
        for i in range(l,9):
            vt.append(0)
        if vt[1] == 0:
            vt[1]=1
        if vt[2] == 0:
            vt[2]=1
        ss=vt[5]-int(vt[5])
        vt[5]=int(vt[5])
        t=time.mktime(tuple(vt))-tz+ss
        vt[5]=vt[5]+ss

    print(vt)
    mjd=v2mjd(vt)
    tt=copy.copy(t)

    if typ == 'UTC':
        t=time.gmtime(t)
    else:
        t=time.localtime(t)

    out=t

    if form == 'string':
        out=time.asctime(t)
    elif form == 'gps':
        out=tt-315964800+leap_seconds(mjd)-9
    elif form == 'unix':
        out=tt
    elif form == 'mjd':
        out=mjd

    return out


def timestruct2v(tstr):
    '''
    time structure to time vector
    '''

    v=[tstr.tm_year,tstr.tm_mon,tstr.tm_mday,tstr.tm_hour,tstr.tm_min,
    tstr.tm_sec,0,0,0]

    return v


def v2mjd(v):
    '''
    vectorial to mjd time
    v   [year,month,day,hour,min,s]
        if month=1, doy can be used instead of day
        s can be float
    '''

    tz=time.timezone
    l=len(v)
    for i in range(l,9):
        v.append(0)
    if v[1] == 0:
        v[1]=1
    if v[2] == 0:
        v[2]=1
    ss=v[5]-int(v[5])
    v[5]=int(v[5])
    
    v=tuple(v)
    t=time.mktime(v)-tz+ss
    # ye=v[1],mo=v[2],da=v[3],ho=v[4],mi=v[5],se=v[6]
    t=t/86400+40587

    return t


def leap_seconds(mjd):
    '''
    leap seconds at mjd
    independent list (see also astropy function)
    we put 0 leapswconds at 1 January 1972
    '''

    leaptimes=[]

    leaptimes.append(41317)
    leaptimes.append(41499)
    leaptimes.append(41683)
    leaptimes.append(42048)
    leaptimes.append(42413)
    leaptimes.append(42778)
    leaptimes.append(43144)
    leaptimes.append(43509)
    leaptimes.append(43874)

    leaptimes.append(44786)
    leaptimes.append(45151)
    leaptimes.append(45516)
    leaptimes.append(46247)
    leaptimes.append(47161)
    leaptimes.append(47892)
    leaptimes.append(48257)
    leaptimes.append(48804)
    leaptimes.append(49169)
    leaptimes.append(49534)
    leaptimes.append(50083)
    leaptimes.append(50630)
    leaptimes.append(51179)
    leaptimes.append(53736)
    leaptimes.append(54832)
    leaptimes.append(56109)
    leaptimes.append(57204)
    leaptimes.append(57754)

    ii=0
    for tt in leaptimes:
        if mjd < tt:
            break
        ii+=1

    return ii




# AstroPy time ------------------------

def set_time(times,form=0):
    '''
    set time values as Time objects
    times   array or list
    form    see https://docs.astropy.org/en/stable/time/index.html#time-format
            if absent it can be interpreted

    The output t are time object, that can be easily converted
    example:
    >>> t = Time('2023-01-01')
    >>> t
    <Time object: scale='utc' format='iso' value=2023-01-01 00:00:00.000>
    >>> t.mjd
    59945.0
    >>> t.gps
    1356566418.0
    '''

    if form == 0:
        t=Time(times)
    else:
        t=Time(times,format=form)

    return t



def t_conv(tin,fmtin,fmtout):
    '''
    time conversion
    tin            input time
    fmtin, fmtout  'str', 'vec', 'mjd', 'gps'
    only fmtout   'obj' time object
                'ut1', 'tt', 'tai',tcb,tcg,tdb
        vec as (2000, 1, 1, 0, 0)
        str as '2000-01-20 00:00:00.000'
    '''

    if fmtin == 'vec':
        aa='datetime.datetime'+str(tin)
        v=eval(aa)
        t=Time(v)
    elif fmtin == 'str':
        t=Time(tin)
    else:
        t=Time(tin,format=fmtin)

    if fmtout == 'mjd':
        t=t.mjd
    elif fmtout == 'gps':
        t=t.gps
    elif fmtout == 'str':
        t=t.iso
    elif fmtout == 'ut1':
        t=t.ut1
    elif fmtout == 'tt':
        t=t.tt
    elif fmtout == 'tai':
        t=t.tai
    elif fmtout == 'tcb':
        t=t.tcb
    elif fmtout == 'tgb':
        t=t.tgb
    elif fmtout == 'tdb':
        t=t.tdb
    elif fmtout == 'vec':
        t=t.value.timetuple()

    return t


def sid_time(tim,loc=0):
    '''
    local or Greenwich (default) sidereal time

    tim   time object or as set_time
    loc   Earth locality object or [lon,lat] (in decimal deg)
    '''

    if not isinstance(tim,Time):
        tim=Time(tim)
    if isinstance(loc,list):
        lon=str(loc[0])+'d'
        lat=str(loc[1])+'d'
        EL=EarLoc(lon,lat)
    elif loc == 0:
        EL=EarLoc(0,0)

    ST=tim.sidereal_time('mean')

    return ST
    


# mjd time ----------------------------

def mjd_now():
    '''
    Return the current GPS time as a float using Astropy.
    '''

    from astropy.time import Time

    return float(Time.now().mjd)



def mjd_phase(t0,long=0,deg=1):
    '''
    phases for physical periods (in deg) or beginning delay
    phisical periods: day, sid, locsid, week
    t0     time
    long   longitude
    deg    = 1  phase in deg
    '''

    sidday=86164.090530833/86400  # at epoch2000

    week=(t0+2)%7
    day=t0-int(t0)
    sid=gmst(t0)*sidday/24
    locsid=gmst(t0,long)*sidday/24

    if deg == 1:
        week=week*360/7
        day=day*360
        sid=sid*360/sidday
        locsid=locsid*360/sidday

    return week,day,sid,locsid





# gps time ----------------------------

def gps_now():
    '''
    Return the current GPS time as a float using Astropy.
    '''

    from astropy.time import Time

    return float(Time.now().gps)


def gps2mjd(tgps): # conversion from gps time to mjd
    pass


# sidereal time ---------------------

def gmst(t,long=0):
    '''
    Greenwich mean sidereal time or local ST
    t     time object or jd or mjd
    long  local longitude in deg or EL object (for local ST)
    '''

    if isinstance(long,EarthLocation):
        aa=long.lon
        long=aa.value
    if isinstance(t,Time):
        t=t.jd
    elif t < 1000000:
        t=t+2400000.5
    
    jd=t
    jd0=int(jd-0.5)+0.5
    h=(jd-jd0)*24

    d=jd-2451545
    d0=jd0-2451545
    T=d/36525

    ST=(6.697374558+0.06570982441908*d0+1.00273790935*h+0.000026*T**2)%24
    ST=ST+long/15

    return ST


def nearest_tsid(ts,t0):
    pass


### ASTRO ###  __________________________________________________

# Constants -----------------------------

def const_table(file=0):
    '''
    file, if it is present, saves the constants table in a file
    '''
    lis0=['G', 'GM_earth', 'GM_jup', 'GM_sun', 'L_bol0', 'L_sun', 'M_earth', 'M_jup', 'M_sun',
     'N_A', 'R', 'R_earth', 'R_jup','R_sun', 'Ryd', 'a0', 'alpha', 'atm', 'au', 'b_wien', 'c', 
     'e', 'eps0', 'g0', 'h', 'hbar','k_B', 'kpc', 'm_e', 'm_n', 'm_p', 'mu0', 'muB', 
     'pc','sigma_T', 'sigma_sb']
    data=[]
    data.append(('Symbol','What','Value','RelErr'))
    for i in lis0:
        ina=eval('const.'+i)
        val=ina.value
        unc=ina.uncertainty
        perc=unc/val
        print(i,ina.name,val,perc)
        data.append((i,ina.name,val,perc))

    st=BASIC.snag_table(data)

    if file != 0:
        BASIC.st2csv(st,file)

    return st


def sites():
    '''
    sites for gravitational antennas
    '''

    virgo_s=EarthLocation.of_site('VIRGO')
    print(virgo_s.lon,virgo_s.lat,virgo_s.height)

    ligol_s=EarthLocation.of_site('LIGO Livingston Observatory')
    print(ligol_s.lon,ligol_s.lat,ligol_s.height)

    ligoh_s=EarthLocation.of_site('LIGO Hanford Observatory')
    print(ligoh_s.lon,ligoh_s.lat,ligoh_s.height)

    kagra_s=EarthLocation.of_site('KAGRA')
    print(kagra_s.lon,kagra_s.lat,kagra_s.height)

    return virgo_s,ligol_s,ligoh_s,kagra_s


# Coordinates --------------------------

def astro_coord(cin, cout, ai, di):
    '''
    ASTRO_COORD   astronomical coordinate conversion, from cin to cout

    cin and cout can be

    'equ'      celestial equatorial: right ascension, declination
    'ecl'      ecliptical: longitude, latitude
    'gal'      galactic longitude, latitude

    ai,di      input coordinates  (deg)
    ao,do      output coordinates (deg)

    epsilon = 23.4392911 deg is the inclination of the Earth orbit on the
    equatorial plane, referred to the standard equinox of year 2000
    (epsilon = 23.4457889 deg, for the standard equinox of year 1950).
    '''

    deg2rad = np.pi/180
    rad2deg = 180/np.pi
    ai = np.array(ai)*deg2rad
    di = np.array(di)*deg2rad

    ao = ai
    do = di

    cai = np.cos(ai)
    sai = np.sin(ai)
    cdi = np.cos(di)
    sdi = np.sin(di)

    # ecliptic
    eps = 23.4392911*deg2rad
    ceps = np.cos(eps)
    seps = np.sin(eps)

    # galactic
    gal = 62.9*deg2rad
    cgal = np.cos(gal)
    sgal = np.sin(gal)
    galnode = 282.85*deg2rad

    if cin == 'hor':
        print('not yet implemented')

    elif cin == 'ecl':
        cdo = np.sqrt((cdi*cai) ** 2+(cdi*sai*ceps-sdi*seps) ** 2)
        sdo = cdi*sai*seps+sdi*ceps
        do = np.arctan(sdo/cdo)
        cao = cdi*cai/np.cos(do)
        sao = (cdi*sai*ceps-sdi*seps)/np.cos(do)
        ao = np.arctan2(sao, cao)

    elif cin == 'gal':
        print('not yet implemented')

    if cout == 'hor':
        print('not yet implemented')

    elif cout == 'ecl':
        cdo = np.sqrt((cdi*cai)**2+(cdi*sai*ceps+sdi*seps)**2)
        sdo = -cdi*sai*seps+sdi*ceps
        do = np.arctan(sdo/cdo)
        cao = cdi*cai/np.cos(do)
        sao = (cdi*sai*ceps+sdi*seps)/np.cos(do)
        ao = np.arctan2(sao, cao)

    elif cout == 'gal':
        print('not yet implemented')

    ao = ao/deg2rad
    do = do/deg2rad

    if ao.shape:
        ao[ao < 0] += 360
    else:
        if ao<0:
            ao += 360

    return ao, do


def astro2rect(a, icrad=0):
    '''
    ASTRO2RECT  conversion from astronomical to rectangular coordinates

    a         astronomical coordinates [ra dec module]
                if a = [ra dec] is converted to [ra dec 1]
    icrad     0 = degrees (default), 1 = radiants
    '''

    if not a[0].shape:
        a = [[a[0]], [a[1]]]

    if len(a) == 2:
        a = a + [[1] * len(a[0])]
       
    a = np.array(a)

    if icrad == 0:
        deg2rad = np.pi/180
        a[0, :] = a[0, :]*deg2rad
        a[1, :] = a[1, :]*deg2rad

    r = np.zeros((3, len(a[0])))
    r[0] = np.cos(a[0, :])*np.cos(a[1, :])*a[2, :]
    r[1] = np.sin(a[0, :])*np.cos(a[1, :])*a[2, :]
    r[2] = np.sin(a[1, :])*a[2, :]
    return r


# AstroPy coordinates ----------------------

def EarLoc(lon,lat,h=0):
    '''
    Earth location object
    lon,lat  angles
    h        height in m

    Angles can be expressed as:

    Angle('10.2345d')               String with 'd' abbreviation for degrees
    Angle('1:2:30.43 degrees')      Sexagesimal degrees
    Angle('1°2′3″')                 Unicode degree, arcmin and arcsec symbols
    Angle('1°2′3″N')                Unicode degree, arcmin, arcsec symbols and direction
    Angle('1d2m3.4s')               Degree, arcmin, arcsec
    if it is expressed as a decimal number xxxx.xx, is equivalent to 'xxxx.xxd'
    '''
 
    if isinstance(lon,(int,float)):
        lon=str(lon)+'d'
    if isinstance(lat,(int,float)):
        lat=str(lat)+'d'

    EL=EarthLocation(lon=lon,lat=lat,height=h)

    return EL


def astropy_coord(cin, cout, ai, di, lon=[], lat=[], time=[]):
    '''
    Astronomical coordinate system change
    Angles are in degrees. Local tsid (in hours) and latitude is needed
    for conversions to and from the horizon coordinates.

    cin and out can be

    'hor'         Horizontal coordinate: azimuth, altitude (az,alt) 
                        YOU SHOULD ADD lat and lon!
    'equ'         celestial equatorial: right ascension, declination (ra,dec)
    'ecl'         ecliptical: longitude, latitude (lon,lat)
    'gal'         galactic longitude, latitude (l,b)

    ai,di         input coordinates
    lon,lat,time  local sidereal time and latitude (only for horizontal coordinates)
                    various form for time see https://docs.astropy.org/en/stable/time/index.html

    Angles can be expressed as:

    Angle('10.2345d')               String with 'd' abbreviation for degrees
    Angle('1:2:30.43 degrees')      Sexagesimal degrees
    Angle('1°2′3″')                 Unicode degree, arcmin and arcsec symbols
    Angle('1°2′3″N')                Unicode degree, arcmin, arcsec symbols and direction
    Angle('1d2m3.4s')               Degree, arcmin, arcsec
    if it is expressed as a decimal number xxxx.xx, is equivalent to 'xxxx.xxd'
    '''

    if isinstance(ai,(int,float)):
        ai=str(ai)+'d'
    if isinstance(di,(int,float)):
        di=str(di)+'d'

    if cin == 'hor':  
        if isinstance(lon,(int,float)):
            lon=str(lon)+'d'
        if isinstance(lat,(int,float)):
            lat=str(lat)+'d'
        time=Time(time)
        EL=EarthLocation(lon=lon,lat=lat)
        obj=AltAz(az=ai,alt=di,location=EL,obstime=time)
    if cin == 'equ':
        obj=SkyCoord(ra=ai,dec=di)
    if cin == 'ecl':
        obj=SkyCoord(l=ai,b=di,frame='geocentrictrueecliptic')
    if cin == 'gal':
        obj=SkyCoord(l=ai,b=di,frame='galactic')

    if cout == 'hor':
        EL=EarthLocation(lon=lon,lat=lat)
        time=Time(time)
        out=obj.transform_to(AltAz(obstime=time,location=EL))
    if cout == 'equ':
        out=obj.icrs
    if cout == 'ecl':
        out=obj.geocentrictrueecliptic
    if cout == 'gal':
        out=obj.galactic

    return out


def ang_sep(lon1,lat1,lon2,lat2):
    '''
    angular separation
    angles in deg
    '''

    lon1=lon1*deg2rad
    lat1=lat1*deg2rad
    lon2=lon2*deg2rad
    lat2=lat2*deg2rad

    out=astropy.coordinates.angular_separation(lon1, lat1, lon2, lat2)

    return out*rad2deg


# Using Skyfield ------------------------

def sf_arrtime(tini,tfin,step):
    '''
    Skyfield time array
    tini    initial time tuple (a,m,d,h,m,s)
    fion    final time tuple
    step    in s
    '''

    ts = load.timescale()

    lini=list(tini)
    l=len(lini)
    for i in range(l,6):
        lini.append(0)
    ini=ts.utc(lini[0],lini[1],lini[2],lini[3],lini[4],lini[5])

    lfin=list(tfin)
    l=len(lfin)
    for i in range(l,6):
        lfin.append(0)
    fin=ts.utc(lfin[0],lfin[1],lfin[2],lfin[3],lfin[4],lfin[5])

    ii=np.arange(0,fin-ini,step/86400)
    lini[2]=lini[2]+ii
    lini[2]=list(lini[2])

    t=ts.utc(lini[0],lini[1],lini[2],lini[3],lini[4],lini[5])

    return t


def sf_coord():
    pass


def sf_earth(t,epht='de440s.bsp'):
    '''
    position and velocity of the Earth
    output:   cartesian position (AU)
            cartesian velocity (in unit of c)
            Einstein effect
    '''

    c=299792458
    G=6.6743e-11
    UA=149597870700
    MS=1.98855e30

    eph = load(epht)
    earth=eph['earth']
    e = earth.at(t)
    p=e.position
    w=e.velocity
    v=w.m_per_s/c
    pp=p.au
    d=pp[1]**2+pp[2]**2+pp[0]**2
    
    rs=2*G*MS/c**2
    einstSol=rs/UA
    einst=einstSol/d

    return pp,v,einst


# Using jpl tables ------------------------

def ant_pos_vel(table_p,table_par,tin,tfi,au=0):
    '''
    position and velocity of an antenna + Einstein effect
        table_p    symbol (e.g.: table_virgo_p)
        table_par  table parameters symbol (e.g.: table_par)
        tin,tfi    vectorial time ini,fin (e.g. [2021,1,20,0])
        au = 0   pos is in light seconds, vel in fraction of c
        au = 1   pos is converted to au
        au = 2   pos is converted to km
    '''

    out,tmjd=GWDATA.extr_doppler(table_p,tin,tfi,table_par)
    pos=out[:,1:4]
    if au == 1:
        pos=pos*299792458/1.495978707e11
    if au == 2:
        pos=pos*299792.458
    vel=out[:,4:7]
    einst=out[:,7]

    return pos,vel,einst,tmjd


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


