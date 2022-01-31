#!/usr/bin/env python
#-*- coding: utf-8 -*-
# encoding: utf-8
from __future__ import unicode_literals

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset
from wrf import to_np, getvar, CoordPair, vertcross
from datetime import timedelta, datetime
import glob, os
import multiprocessing
import time
from scipy.interpolate import griddata
from  matplotlib import colors, patches
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import LinearSegmentedColormap

def main(xpto, ypto, lab):

    print(lab)

    # CPU cores
    PROCESS_LIMIT = 8

    os.system('rm -f out*.txt')

    lista = glob.glob("wrfout_d03*:00")
    lista.sort()

    nc_test = Dataset(lista[0], "r")
    points = np.array( (nc_test.variables['XLONG'][0].flatten(), nc_test.variables['XLAT'][0].flatten()) ).T

    # Height of the station to calculate MSLP
    hgt_example = np.array(griddata( points, nc_test.variables['HGT'][0].flatten(), (np.array([xpto]), np.array([ypto])) ))
    nc_test.close()

    # interpoling data
    for wrf in lista:

        finterp(wrf, points, xpto, ypto, hgt_example)
#        process=multiprocessing.Process(target=finterp,args=(wrf, points, xpto, ypto, hgt_example))
#        while(len(multiprocessing.active_children()) == PROCESS_LIMIT):
#            time.sleep(1)
#        process.start()

    # colecting data
    # t2m, mslp, rh2m,wspeed, wind_dir_cardinal, u_wind_msI, v_wind_msI, TEE, nubosidad, alturas, rain

    itemp2m = []
    imslp   = []
    irh2m   = []
    iwspeed = []
    iwdirec = []
    iuwind  = []
    ivwind  = []
    itee    = []
    iclouds = []
    iheigths = []
    irain   = []

    times = []

    c=0
    for wrf in lista:

        times.append(parse_date(wrf.split('_')[-2]+' '+wrf.split('_')[-1]))
        tmpdata = np.loadtxt('out_'+wrf+'.txt').flatten()

        itemp2m.append(tmpdata[0])
        imslp.append(tmpdata[1])
        irh2m.append(tmpdata[2])
        iwspeed.append(tmpdata[3])
        iwdirec.append(tmpdata[4])
        iuwind.append(tmpdata[5])
        ivwind.append(tmpdata[6])
        itee.append(tmpdata[7])

        if c != 0:
            tmpdata0 = np.loadtxt('out_'+lista[c-1]+'.txt').flatten()
            tmpdata1 = np.loadtxt('out_'+lista[c]+'.txt').flatten()
            irain.append(np.abs(tmpdata1[8]-tmpdata0[8]))
        else:
            irain.append(0.0)

        tmpdata2 = np.loadtxt('out2_'+wrf+'.txt')

        iclouds.append(tmpdata2[:,0])
        iheigths.append(tmpdata2[:,1])

        c += 1


    # set the starttime and endtime for plotting, 24 hour range
    endtime = parse_date(lista[-1].split('_')[-2]+' '+lista[-1].split('_')[-1])
    starttime = parse_date(lista[0].split('_')[-2]+' '+lista[0].split('_')[-1])

    endtime1 = parse_date2(lista[-1].split('_')[-2]+' '+lista[-1].split('_')[-1])
    starttime1 = parse_date2(lista[0].split('_')[-2]+' '+lista[0].split('_')[-1])

    totalhours = diffdates(lista[0].split('_')[-2]+' '+lista[0].split('_')[-1],lista[-1].split('_')[-2]+' '+lista[-1].split('_')[-1])
    dthours    = diffdates(lista[0].split('_')[-2]+' '+lista[0].split('_')[-1],lista[1].split('_')[-2]+' '+lista[1].split('_')[-1])


    #-------------------------------------------------------------------------------------------------

    f, axarr = plt.subplots(6, figsize=(9, 18), dpi=300, sharex=False)
    f.subplots_adjust(hspace=0.5)

    #-------------------------------------------------------------------------------------------------

    mn1, mx1 = (np.min(itemp2m),np.max(itemp2m))
    mn2, mx2 = (np.min(itee),np.max(itee))
    mn0, mx0 = (np.min(np.array([mn1, mn2])), np.max(np.array([mx1, mx2])))

    cmap1,clevs1=cm_temp()
    X,Y = np.meshgrid(np.arange(0,len(times),1), clevs1)
    axarr[0].contourf(X, Y, Y, clevs1, cmap=cmap1,zorder=2)
    axarr[0].fill_between(np.arange(0,len(times),1), itemp2m, mx0, alpha=1, facecolor='w', zorder=3)

    axarr[0].set_title('Pronóstico Numérico del modelo WRF-SisPI sobre el comportamiento de las variables meteorológicas \npara las próximas '+totalhours+' horas en el municipio '+lab.upper()+'\nInicializado '+starttime1+' y Válido hasta '+endtime1+'\nTemperatura a 2 metros [℃] / Temp. Efectiva Equivalente (Sensación Térmica)'.upper(), fontsize = 10, fontweight="bold")
    # '\nCoordenadas geog.: '+str(np.round(ypto,4))+'N/'+str(np.round(xpto,4))+'W\n'

    ln1 = axarr[0].plot(np.arange(0,len(times),1), itemp2m, '-', color='k',linewidth=2, markeredgewidth=2, markersize=2, alpha=0.6, label='Temp.')
#    ln1 = axarr[0].plot_date(times, itemp2m, '-', color='k',linewidth=1, markeredgewidth=0, markersize=0, alpha=0.6, label='Temp.')
#    axarr[0].fill_between(times, itemp2m, 0, alpha=0.1, facecolor='b')
#    axarr[0].set_xlim(starttime, endtime)
    axarr[0].set_xlim(np.arange(0,len(times),1)[0], np.arange(0,len(times),1)[-1])
    axarr[0].set_ylim(mn0, mx0)
#    axarr[0,0].grid()#(b=True, which='major', axis='y', color='k', linestyle='--', linewidth=0.5)
    axarr[0].set_ylabel('Temperatura [℃]', multialignment='center', fontsize = 9, color='blue')
    axarr[0].set_xticks(np.arange(0,len(times),1)[::2])
    axarr[0].set_xticklabels( times[::2], rotation=0)


    ax0 = axarr[0].twinx()
    ln2 = ax0.plot(np.arange(0,len(times),1), itee, 'o', color='k',linewidth=0, markersize=4.0, markerfacecolor='w', markeredgecolor='k', markeredgewidth=1.0, label='Sens.Term.')
#    ln2 = ax0.plot_date(times, itee, 'o', color='k',linewidth=0, markeredgewidth=2, markersize=2, label='Sens.T.')
    ax0.set_ylabel('Sens. Térmica', multialignment='center', fontsize = 9, color='blue')
    ax0.set_xlim(np.arange(0,len(times),1)[0], np.arange(0,len(times),1)[-1])
    ax0.set_ylim(mn0, mx0)
    ax0.set_xticks(np.arange(0,len(times),1)[::2])
    ax0.set_xticklabels( times[::2], rotation=0)
    ax0.grid(which='major', color='#CCCCCC', linestyle='--')

    ax0.legend(loc='lower left', fontsize=7)

    #-------------------------------------------------------------------------------------------------

    mn2, mx2 = (np.min(irh2m),np.max(irh2m))

    cmap1,clevs1=cm_rh2m()
    X,Y = np.meshgrid(np.arange(0,len(times),1), clevs1)
    axarr[1].contourf(X, Y, Y, clevs1, cmap=cmap1,zorder=2)
    axarr[1].fill_between(np.arange(0,len(times),1), irh2m, mx2, alpha=1, facecolor='w', zorder=3)

    axarr[1].set_title('Humedad Relativa [%]'.upper(), fontsize = 10, fontweight="bold")

    ln1 = axarr[1].plot(np.arange(0,len(times),1), irh2m, '-', color='k',linewidth=1, markeredgewidth=0, markersize=0, alpha=0.6, label='Temp.')
#    ln1 = axarr[0].plot_date(times, irh2m, '-', color='k',linewidth=1, markeredgewidth=0, markersize=0, alpha=0.6, label='Temp.')
#    axarr[0].fill_between(times, irh2m, 0, alpha=0.1, facecolor='b')
#    axarr[0].set_xlim(starttime, endtime)
    axarr[1].set_xlim(np.arange(0,len(times),1)[0], np.arange(0,len(times),1)[-1])
    axarr[1].set_ylim(mn2, mx2)
    axarr[1].grid()#(b=True, which='major', axis='y', color='k', linestyle='--', linewidth=0.5)
    axarr[1].set_ylabel('Humedad Rel. [%]', multialignment='center', fontsize = 9, color='blue')
    axarr[1].set_xticks(np.arange(0,len(times),1)[::2])
    axarr[1].set_xticklabels( times[::2], rotation=0)

    ax1 = axarr[1].twinx()
    ln2 = ax1.plot(np.arange(0,len(times),1), irh2m, '-', color='k',linewidth=1, markeredgewidth=2, markersize=2, label='Hum. Rel.')
#    ln2 = ax1.plot_date(times, imslp, 'o', color='k',linewidth=0, markeredgewidth=2, markersize=2, label='Sens.T.')
    ax1.set_ylabel('Humedad Rel. [%]', multialignment='center', fontsize = 9, color='blue')
    ax1.set_xlim(np.arange(0,len(times),1)[0], np.arange(0,len(times),1)[-1])
    ax1.set_ylim(mn2, mx2)
    ax1.set_xticks(np.arange(0,len(times),1)[::2])
    ax1.set_xticklabels( times[::2], rotation=0)
    ax1.grid(which='major', color='#CCCCCC', linestyle='--')

#    ax1.legend(loc='lower left', fontsize=6)

    #-------------------------------------------------------------------------------------------------

    mn1, mx1 = (np.min(imslp),np.max(imslp))

    cmap1,clevs1=cm_mslp()
    X,Y = np.meshgrid(np.arange(0,len(times),1), clevs1)
    axarr[2].contourf(X, Y, Y, clevs1, cmap=cmap1,zorder=2)
    axarr[2].fill_between(np.arange(0,len(times),1), imslp, mx1, alpha=1, facecolor='w', zorder=3)

    axarr[2].set_title('Presión Ajustada el Nivel Medio del Mar [hPa]'.upper(), fontsize = 10, fontweight="bold")
    ln1 = axarr[2].plot(np.arange(0,len(times),1), imslp, '-', color='k',linewidth=1, markeredgewidth=0, markersize=0, alpha=0.6, label='Temp.')
#    ln1 = axarr[2].plot_date(times, iwspeed, '-', color='k',linewidth=1, markeredgewidth=0, markersize=0, alpha=0.6, label='Temp.')
#    axarr[2].fill_between(times, iwspeed, 0, alpha=0.1, facecolor='b')
#    axarr[2].set_xlim(starttime, endtime)
    axarr[2].set_xlim(np.arange(0,len(times),1)[0], np.arange(0,len(times),1)[-1])
    axarr[2].set_ylim(mn1, mx1)
    axarr[2].grid()#(b=True, which='major', axis='y', color='k', linestyle='--', linewidth=0.5)
    axarr[2].set_ylabel('Presión Atm. [hPa]', multialignment='center', fontsize = 9, color='blue')
    axarr[2].set_xticks(np.arange(0,len(times),1)[::2])
    axarr[2].set_xticklabels( times[::2], rotation=0)

    ax2 = axarr[2].twinx()
    ln2 = ax2.plot(np.arange(0,len(times),1), imslp, '-', color='k',linewidth=1, markeredgewidth=2, markersize=2, label='Presión Atm.')
#    ln2 = ax2.plot_date(times, imslp, 'o', color='k',linewidth=0, markeredgewidth=2, markersize=2, label='Sens.T.')
    ax2.set_ylabel('Presión Atm. [hPa]', multialignment='center', fontsize = 9, color='blue')
    ax2.set_xlim(np.arange(0,len(times),1)[0], np.arange(0,len(times),1)[-1])
    ax2.set_ylim(mn1, mx1)
    ax2.set_xticks(np.arange(0,len(times),1)[::2])
    ax2.set_xticklabels( times[::2], rotation=0)
    ax2.grid(which='major', color='#CCCCCC', linestyle='--')

#    ax2.legend(loc='lower left', fontsize=6)

    #-------------------------------------------------------------------------------------------------

    mn1, mx1 = (np.min(iwspeed),np.max(iwspeed))

    cmap1,clevs1=cm_wind()
    X,Y = np.meshgrid(np.arange(0,len(times),1), clevs1)
    axarr[3].contourf(X, Y, Y, clevs1, cmap=cmap1,zorder=2)
    axarr[3].fill_between(np.arange(0,len(times),1), iwspeed, mx1, alpha=1, facecolor='w', zorder=3)

    axarr[3].set_title('Velocidad y direccón del viento [km/h]'.upper(), fontsize = 10, fontweight="bold")
    ln1 = axarr[3].plot(np.arange(0,len(times),1), iwspeed, '-', color='k',linewidth=1, markeredgewidth=0, markersize=0, alpha=0.6, label='Temp.')
#    ln1 = axarr[0,1].plot_date(times, iwspeed, '-', color='k',linewidth=1, markeredgewidth=0, markersize=0, alpha=0.6, label='Temp.')
#    axarr[0,1].fill_between(times, iwspeed, 0, alpha=0.1, facecolor='b')
#    axarr[0,1].set_xlim(starttime, endtime)
    axarr[3].set_xlim(np.arange(0,len(times),1)[0], np.arange(0,len(times),1)[-1])
    axarr[3].set_ylim(mn1, mx1)
    axarr[3].grid()#(b=True, which='major', axis='y', color='k', linestyle='--', linewidth=0.5)
    axarr[3].set_ylabel('Vel. Viento. [km/h]', multialignment='center', fontsize = 9, color='blue')
    axarr[3].set_xticks(np.arange(0,len(times),1)[::2])
    axarr[3].set_xticklabels( times[::2], rotation=0)

    ywind = np.arange(0,len(times),1)
    ywind[:] = (mx1 - mn1)/2

    ax3 = axarr[3].twinx()
    ln2 = ax3.plot(np.arange(0,len(times),1), iwspeed, '-', color='k',linewidth=1, markeredgewidth=2, markersize=2, label='Presión Atm.')
#    ln2 = ax3.plot_date(times, imslp, 'o', color='k',linewidth=0, markeredgewidth=2, markersize=2, label='Sens.T.')
    ax3.quiver(np.arange(0,len(times),1), ywind, np.array(iuwind)/iwspeed, np.array(ivwind)/iwspeed, alpha=0.7, scale_units='inches', scale=1.5)#, minlength=0)

    ax3.set_ylabel('Vel. Viento. [km/h]', multialignment='center', fontsize = 9, color='blue')
    ax3.set_xlim(np.arange(0,len(times),1)[0], np.arange(0,len(times),1)[-1])
    ax3.set_ylim(mn1, mx1)
    ax3.set_xticks(np.arange(0,len(times),1)[::2])
    ax3.set_xticklabels( times[::2], rotation=0)
    ax3.grid(which='major', color='#CCCCCC', linestyle='--')

#    ax3.legend(loc='lower left', fontsize=6)

    #-------------------------------------------------------------------------------------------------

    mn1, mx0 = (np.min(irain),np.max(irain))
    mx1 = mx0 + np.abs(mx0-mn1)/4

    if dthours == '1':
        tprecip='1 hora'
    else:
        tprecip=dthours+' horas'

    axarr[4].set_title('Precipitación [mm/'.upper()+tprecip+']'.upper(), fontsize = 10, fontweight="bold")
    ln1 = axarr[4].bar(np.arange(0,len(times),1), irain, color='blue', alpha=0.6, label='Precip.')
    axarr[4].set_xticks(np.arange(0,len(times),1)[::2])
    axarr[4].set_xticklabels( times[::2], rotation=0)

    autolabel(ln1, axarr[4])

    axarr[4].set_xlim(np.arange(0,len(times),1)[0], np.arange(0,len(times),1)[-1])
    axarr[4].set_ylim(mn1, mx1)
    axarr[4].grid()#(b=True, which='major', axis='y', color='k', linestyle='--', linewidth=0.5)
    axarr[4].set_ylabel('Precipitación [mm/'+tprecip+']', multialignment='center', fontsize = 9, color='blue')

    ax4 = axarr[4].twinx()
    ax4.set_ylabel('Precipitación [mm/'+tprecip+']', multialignment='center', fontsize = 9, color='blue')
    ax4.set_xlim(np.arange(0,len(times),1)[0], np.arange(0,len(times),1)[-1])
    ax4.set_ylim(mn1, mx1)
    ax4.set_xticks(np.arange(0,len(times),1)[::2])
    ax4.set_xticklabels( times[::2], rotation=0)


    #-------------------------------------------------------------------------------------------------

    mn1, mx1 = (np.min(iheigths)/1000,np.max(iheigths)/1000)

    cmap1,clevs1=cm_cloud()
    cloudarray, heigths, X = (np.zeros(shape=(len(iheigths[0]), len(np.arange(0,len(times),1)))), np.zeros(shape=(len(iheigths[0]), len(np.arange(0,len(times),1)))), np.zeros(shape=(len(iheigths[0]), len(np.arange(0,len(times),1)))))
    for i in range(len(np.arange(0,len(times),1))):
        cloudarray[:,i] = iclouds[i][:]
        heigths[:,i]    = iheigths[i][:]
        X[:,i]          = i

#    print X, heigths, cloudarray, np.max(cloudarray)
    CLFRA = axarr[5].contourf(X, heigths/1000, cloudarray*100, clevs1, cmap=cmap1, alpha=0.99, extend="max")  # , zorder=2
    CLFRA.cmap.set_under((0.0, 0.0, 0.0))
    CLFRA.cmap.set_over((1.0, 1.0, 1.0))

    axarr[5].set_title('Perfíl Vertical de Nubosidad [%]'.upper(), fontsize = 10, fontweight="bold")

    axarr[5].set_xlim(np.arange(0,len(times),1)[0], np.arange(0,len(times),1)[-1])
    axarr[5].set_ylim(mn1, mx1)

    axarr[5].grid()#(b=True, which='major', axis='y', color='k', linestyle='--', linewidth=0.5)
    axarr[5].set_ylabel('Altura sobre el nivel del suelo [km]', multialignment='center', fontsize = 9, color='blue')

    axarr[5].set_xticks(np.arange(0,len(times),1)[::2])
    axarr[5].set_xticklabels( times[::2], rotation=0)


    ax5 = axarr[5].twinx()
    ax5.set_ylabel('Altura sobre el nivel del suelo [km]', multialignment='center', fontsize = 9, color='blue')
    ax5.set_xlim(np.arange(0,len(times),1)[0], np.arange(0,len(times),1)[-1])
    ax5.set_ylim(mn1, mx1)
    ax5.set_xticks(np.arange(0,len(times),1)[::2])
    ax5.set_xticklabels( times[::2], rotation=0)


    #-------------------------------------------------------------------------------------------------
#    f.suptitle('Comportamiento de las variables meteorológicas para las próximas '.upper()+totalhours+' horas en el municipio '.upper()+lab.upper()+' \n Pronóstico Numérico del modelo WRF-SisPI Inicializado '.upper()+starttime1+' y Válido hasta '.upper()+endtime1+'\nCoordenadas geográficas: '.upper()+str(np.round(ypto,4))+'N / '+str(np.round(xpto,4))+'W', fontsize = 10)

    plt.savefig('meteogram_'+lab+'.png', dpi=300, bbox_inches='tight', pad_inches=0)
    plt.close('all')


def autolabel(rects, ax):
    # attach some text labels
    for rect in rects:
        height = rect.get_height()
#        print height
        if height > 0.2:
#            ax.text(rect.get_x()+rect.get_width()/2., 1.05*height, '%d'%int(height), ha='center', va='bottom', fontsize=8)
            ax.text(rect.get_x()+rect.get_width()/2., 1.05*height, np.round(height,1), ha='center', va='bottom', fontsize=8)


def diffdates(date1, date2):
#    utctime1 = datetime.strptime(date1.decode('ascii'), '%Y-%m-%d %H:%M:%S')
#    utctime2 = datetime.strptime(date2.decode('ascii'), '%Y-%m-%d %H:%M:%S')
    utctime1 = datetime.strptime(date1, '%Y-%m-%d %H:%M:%S')
    utctime2 = datetime.strptime(date2, '%Y-%m-%d %H:%M:%S')

    difft = abs(utctime2 - utctime1).total_seconds()/3600

    return str(int(difft))


def parse_date(date):
#    utctime = datetime.strptime(date.decode('ascii'), '%Y-%m-%d %H:%M:%S')
    utctime = datetime.strptime(date, '%Y-%m-%d %H:%M:%S')

    if str(utctime)[5:7] in ['11','12','01','02','03','04']:
        utcdif = 5
    if str(utctime)[5:7] in ['05','06','07','08','09','10']:
        utcdif = 4

    localtime = utctime - timedelta(hours=utcdif)

#    print str(localtime.strftime("%d %I:%M %p"))
    return str(localtime.strftime("%I%p\n%d"))

def parse_date2(date):
#    utctime = datetime.strptime(date.decode('ascii'), '%Y-%m-%d %H:%M:%S')
    utctime = datetime.strptime(date, '%Y-%m-%d %H:%M:%S')

    if str(utctime)[5:7] in ['11','12','01','02','03','04']:
        utcdif = 5
    if str(utctime)[5:7] in ['05','06','07','08','09','10']:
        utcdif = 4

    localtime = utctime - timedelta(hours=utcdif)

    return str(localtime.strftime("%Y-%m-%d %I:%M %p"))


def finterp(wrf, points, xpto, ypto, hgt_example):

#    if not os.path.exists('out2_'+wrf+'.txt'):

        nc = Dataset(wrf, "r")

        # Rain

        rainArrayC =  nc.variables['RAINC'][0]
        rainArrayNC =  nc.variables['RAINNC'][0]
        rainArray = rainArrayC + rainArrayNC

        rain = griddata( points, np.abs(rainArray).flatten(), (np.array([xpto]), np.array([ypto])) )

        # Cloud fraction

        z = getvar(nc, "z")
        cloudfr =  getvar(nc, "CLDFRA")
        hlevels =  getvar(nc, "height")

        start_point = CoordPair(lat=float(ypto), lon=float(xpto))
        end_point = CoordPair(lat=float(ypto), lon=float(xpto)+0.1)
        cloudfr_cross = vertcross(cloudfr, z, wrfin=nc, start_point=start_point, end_point=end_point, latlon=True, meta=True)
        hlevels_cross = vertcross(hlevels, z, wrfin=nc, start_point=start_point, end_point=end_point, latlon=True, meta=True)

        nubosidad, alturas  = (to_np(cloudfr_cross)[1:,0], to_np(hlevels_cross)[1:,0])

        # Grab these variables for now
        tempsI = griddata( points, nc.variables['T2'][0].flatten(), (np.array([xpto]), np.array([ypto])) )
        psfcI = griddata( points, nc.variables['PSFC'][0].flatten(), (np.array([xpto]), np.array([ypto])) )
        qhumI = griddata( points, nc.variables['Q2'][0].flatten(), (np.array([xpto]), np.array([ypto])) )

        # Mean Sea Level Presion
        stemps = tempsI+6.5*hgt_example/1000.
        mslp = psfcI*np.exp(9.81/(287.0*stemps)*hgt_example)*0.01 + (6.7 * hgt_example / 1000)

        # Dew point Temperature
        es = 6.112 * np.exp(17.67 * tempsI/(tempsI + 243.5))
        w = qhumI/(1-qhumI)
        e = (w * psfcI / (.622 + w)) / 100
        td2m = (243.5 * np.log(e/6.112))/(17.67-np.log(e/6.112))

        # Air Temperature
        t2m = tempsI - 273.15

        # Relative Humidity
        ens=6.1*10*(7.5*td2m/(237.7+td2m))
        esa=6.1*10*(7.5*t2m/(237.7+t2m))
        rh2m = 100. * (ens/esa)

        # Wind Speed
        u_wind_msI = griddata( points, nc.variables['U10'][0].flatten(), (np.array([xpto]), np.array([ypto])) )
        v_wind_msI = griddata( points, nc.variables['V10'][0].flatten(), (np.array([xpto]), np.array([ypto])) )
        wind_abs = np.sqrt(u_wind_msI*u_wind_msI + v_wind_msI*v_wind_msI)
        wspeed = wind_abs*3.6

        # Wind Direction
#        wind_dir_trig_to = np.arctan2(u_wind_msI/wind_abs, v_wind_msI/wind_abs) 
#        wind_dir_trig_to_degrees = wind_dir_trig_to * 180/np.pi
#        wind_dir_trig_from_degrees = wind_dir_trig_to_degrees + 180
#        wind_dir_cardinal = 90 - wind_dir_trig_from_degrees

        wind_dir_cardinal =  griddata( points, to_np(getvar(nc, "wspd_wdir10"))[1,:,:].flatten(), (np.array([xpto]), np.array([ypto])))

        # Temperatura Efectiva Equivalente
        difT = t2m-37
        TE=t2m-(100-rh2m)/80*(0.00439*(difT*difT)+0.456*difT+9.5)
        TEE=TE+(wind_abs*0.67)*((0.11*difT-0.13)-0.002*difT*(100-rh2m))


        # Salida a ficheros txt
        outlist = [t2m, mslp, rh2m,wspeed, wind_dir_cardinal, u_wind_msI, v_wind_msI, TEE, rain]
        outlist = np.array(outlist)
        np.savetxt('out_'+wrf+'.txt',outlist)

        out2 = np.zeros(shape=(len(nubosidad),2))
        out2[:,0] = nubosidad
        out2[:,1] = alturas
        np.savetxt('out2_'+wrf+'.txt',out2)


        nc.close()

def cm_cloud():

#    clevs1 = np.linspace(0,100,25)
#    cmap1 = plt.cm.gray_r

#    clevs1 = np.array([0,5,15,25,50,75,90,100])
#    windcolors = ((255,255,255),(238,238,238),(225,225,225),(212,212,212),(200,200,200),(170,170,170),(130,130,130),(90,90,90))

    clevs1 = np.array([0,5,20,40,65,85,95,100])
    windcolors = ((255,255,255),(210,210,210),(190,190,190),(170,170,170),(150,150,150),(125,125,125),(100,100,100),(50,50,50))
    cmap1 = LinearSegmentedColormap.from_list("wind", cnv_to_rgb(windcolors), N=len(windcolors), gamma=1.0)

    return cmap1, clevs1

def cm_temp():
    a = np.array([0,2,4,6,8,10,12,14,16,17,18,19,19.5,20,20.5,21.0,21.5,22.0,22.5,23.0,23.5,24.0,25.0,25.5,26.0,26.5,27.0,27.5,28.0,28.5,29.0,29.7,30.7,31.5,32,33,34,35,36])

    # Bins normalized between 0 and 1
    norm = [(float(i)-min(a))/(max(a)-min(a)) for i in a]

    C = np.array([[140,140,140],
        [188,188,188],
        [230,230,230],
        [255,255,255],
        [190,190,255],
        [160,140,255],
        [112,96,220],
        [90,70,200],
        [55,40,165],
        [20,0,130],
        [20,100,210],
        [40,130,240],
        [80,165,245],
        [10,145,70],
        [40,170,85],
        [75,175,105],
        [120,195,110],
        [150,205,125],
        [190,225,135],
        [200,230,140],
        [220,240,150],
        [240,240,195],
        [240,235,140],
        [240,215,130],
        [245,200,90],
        [240,175,75],
        [230,155,60],
        [240,135,45],
        [225,115,0],
        [250,80,60],
        [240,15,105],
        [140,0,0],
        [190,0,0],
        [100,0,5],
        [120,80,70],
        [140,100,90],
        [180,140,130],
        [225,190,180],
        [248,219,214]])/255.

    # Create a tuple for every color indicating the normalized position on the colormap and the assigned color.
    COLORS = []
    for i, n in enumerate(norm):
        COLORS.append((n, C[i]))

    cmap = colors.LinearSegmentedColormap.from_list("Temperature", COLORS)

    return cmap,a


def cm_rh2m():
    a = np.array([0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100])

    # Bins normalized between 0 and 1
    norm = [(float(i)-min(a))/(max(a)-min(a)) for i in a]

    C = np.array([[245,245,245],
        [255,255,160],
        [254,254,99],
        [244,244,110],
        [255,210,35],
        [255,163,25],
        [255,89,25],
        [230,122,101],
        [237,145,124],
        [239,178,146],
        [248,199,178],
        [255,230,230],
        [215,225,255],
        [150,210,255],
        [50,190,255],
        [20,150,255],
        [9,118,240],
        [4,75,185],
        [0,40,140],
        [0,128,0],
        [59,179,59],
        [118,230,118]])/255.

    # Create a tuple for every color indicating the normalized position on the colormap and the assigned color.
    COLORS = []
    for i, n in enumerate(norm):
        COLORS.append((n, C[i]))

    cmap = colors.LinearSegmentedColormap.from_list("Temperature", COLORS)

    return cmap,a

def cm_mslp():

    a = [992,994,996,998,1000,1002,1004,1006,1007,1008,1009,1010,1011,1012,1013,1014,1015,1016,1017,1018,1019,1020,1022,1024,1026,1028,1030,1032]
#    a = [800,850,900,950,970,990,1000,1002,1004,1006,1008,1010,1012,1013,1014,1015,1016,1017,1018,1019,1020,1022,1025,1028,1030,1035,1040,1050]
    clevs = np.array(a)

    # Normalize the bin between 0 and 1 (uneven bins are important here)
    norm = [(float(i)-min(a))/(max(a)-min(a)) for i in a]

    C = np.array([[100,0,0],
                    [170,0,0],
                    [215,0,0],
                    [255,0,0],
                    [255,30,30],
                    [255,70,70],
                    [255,90,90],
                    [255,110,110],
                    [255,135,135],
                    [255,175,175],
                    [255,190,190],
                    [255,210,210],
                    [255,230,230],
                    [255,240,240],
                    [255,255,255],
                    [240,240,255],
                    [230,230,255],
                    [210,210,255],
                    [190,190,255],
                    [175,175,255],
                    [135,135,255],
                    [110,110,255],
                    [90,90,255],
                    [70,70,255],
                    [30,30,255],
                    [0,0,255],
                    [0,0,215],
                    [0,0,170],
                    [0,0,100]])

    # Create a tuple for every color indicating the normalized position on the colormap and the assigned color.
    COLORS = []
    for i, n in enumerate(norm):
        COLORS.append((n, np.array(C[i])/255.))

    # Create the colormap
    cmap = colors.LinearSegmentedColormap.from_list("mslp", COLORS)

    return cmap,clevs


def cm_wind():

    windcolors = ((235,235,235),(215,225,255),(181,201,255),(142,178,255),(127,150,255),(99,112,248),(0,99,255),(0,100,210),
    (0,150,150),(0,160,70),(0,198,51),(50,225,25),(99,235,0),(140,255,0),(198,255,51),(230,255,0),(255,245,0),(255,220,0),
    (255,188,0),(255,125,0),(255,85,0),(255,0,0),(215,0,0),(170,0,0),(105,0,70),(170,0,100),(240,0,130),(240,0,160),
    (245,120,190),(250,190,230),(255,230,235),(255,251,253))
    clevs1 = np.array([0.0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,115,120,125,130,135,140,145,150])
    cmap1 = LinearSegmentedColormap.from_list("wind", cnv_to_rgb(windcolors), N=len(windcolors), gamma=1.0)

    return cmap1,clevs1


def cnv_to_rgb(clist):

    newcolors = []
    for i in range(len(clist)):
        newcolors.append((float(clist[i][0])/255,float(clist[i][1])/255,float(clist[i][2])/255))

    return newcolors




if __name__ == '__main__':

#    Estaciones insmet

#    stn = [[310,-84.95144722222223,21.86636944444445],
#    [312,-83.95616611509935,22.65925530718154],
#    [313,-84.11447266378347,22.14835263558624],
#    [314,-83.82962097303847,22.27795673844591],
#    [315,-83.65362274962783,22.40472916174542],
#    [317,-83.30684041532359,22.56350373334215],
#    [316,-83.54447499999999,22.76838333333333],
#    [318,-83.17009166666666,22.932625],
#    [320,-82.51251111111111,22.77972777777778],
#    [376,-82.53784166666667,22.98218611111111],
#    [322,-82.29378055555556,22.73920277777778],
#    [323,-82.03909444444444,22.85173611111111],
#    [340,-81.93792500000001,23.00235555555556],
#    [374,-82.14043611111111,22.99865555555556],
#    [375,-82.11507777777777,22.78171388888889],
#    [325,-82.341275,23.14361388888889],
#    [373,-82.38856666666668,22.97857777777778],
#    [327,-81.54416666666667,22.76714722222222],
#    [328,-81.23101944444444,23.15935555555555],
#    [329,-81.01693333333333,22.81716111111111],
#    [330,-81.14395561133796,22.79717337190929],
#    [331,-81.13865000000001,22.53484722222223],
#    [332,-80.92548333333333,22.68476111111111],
#    [333,-81.03238888888889,22.07131111111111],
#    [344,-80.44453611111112,22.18983888888889],
#    [335,-80.82605,22.37320833333333],
#    [308,-79.98054722222223,22.11431944444445],
#    [326,-80.22605277777778,22.58687777777778],
#    [338,-80.09233055555555,22.80502222222222],
#    [343,-79.99189444444444,22.46256944444444],
#    [348,-79.47074166666667,22.49694722222222],
#    [341,-79.23263888888889,21.73954722222223],
#    [342,-80.01520277777777,21.91815555555556],
#    [349,-79.44831666666667,21.96981388888889],
#    [337,-79.98935277777778,21.78257222222222],
#    [339,-78.36786944444444,22.53857222222223],
#    [345,-78.85530555555555,21.62448888888889],
#    [346,-78.79763888888888,21.76073055555555],
#    [347,-78.77055833333333,22.1634],
#    [350,-78.25039166666667,21.51733333333333],
#    [351,-78.00082777777777,20.7375],
#    [352,-78.11909166666666,21.84131111111111],
#    [353,-77.24747499999999,21.55975277777778],
#    [354,-77.32135,21.14673611111111],
#    [355,-77.84985833333333,21.42277222222223],
#    [357,-76.94372222222222,20.94478333333333],
#    [358,-76.61633055555555,21.20903055555555],
#    [362,-76.5403361111111,20.9344],
#    [365,-75.62074166666666,21.07167222222222],
#    [370,-75.78237222222222,20.67265555555556],
#    [371,-75.79000555555555,20.48724166666667],
#    [372,-76.22098611111112,20.88517222222222],
#    [378,-76.34462777777777,21.06247222222222],
#    [359,-77.16124444444445,20.31126111111111],
#    [360,-77.72213055555555,19.84056666666666],
#    [361,-76.90107500000001,20.69467777777778],
#    [377,-76.89650555555556,20.34693888888889],
#    [363,-76.26194722222222,20.29731666666667],
#    [364,-75.81717499999999,20.04443888888889],
#    [319,-74.84407499999999,20.14784722222222],
#    [334,-74.95652222222222,20.36753888888889],
#    [356,-74.44640277777778,20.2979],
#    [368,-75.23387777777778,20.1348],
#    [369,-74.14604745735608,20.24167470608382],
#    [309,-82.84908333333333,21.83595555555555],
#    [321,-82.75730833333333,21.73333055555555],
#    [366,-75.62634722222222,20.011625],
#    [324,-82.55784444444444,21.5626]]


#    # Puntos de interes de occidente
    xptos = [-82.3858,-82.1519,-82.7614,-83.6987,-81.5787,-81.1280,-79.9662,-80.4526,-82.5044]
    yptos = [ 23.1245, 22.9615, 22.8107, 22.4167, 23.0477, 22.5260, 22.4071, 22.1457, 22.8907]
    xlab  = ['Plaza_(La_Habana)', 'San_Jose_(Mayabeque)', 'Artemisa_(Artemisa)', 'Pinar_del_Rio_(Pinar_del_R.)','Matanzas_(Matanzas)','Jaguey_Grande_(Matanzas)','Santa_Clara_(Villa_Clara)','Cienfuegos_(Cienfuegos)','San_Antonio_Banos_(Artemisa)']

#    for xpto, ypto, lab in zip(xptos,yptos,xlab):

#        main(xpto, ypto, lab)

#    xptos = [-80.5, -81.5]
#    yptos = [23., 22.0]
#    xlab  = ['a','b']

    for xpto, ypto, lab in zip(xptos,yptos,xlab):

        main(xpto, ypto, lab)


