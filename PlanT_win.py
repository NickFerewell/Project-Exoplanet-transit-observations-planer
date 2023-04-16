import numpy as np
import matplotlib as mpl
##mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import gridspec
import os
import shutil

from string import Template
import urllib.parse
import ephem

from astropy.time import Time, TimeDelta
from astropy.io import ascii
import astropy.units as u
from astropy.coordinates import get_sun, get_moon, EarthLocation, AltAz, SkyCoord
from astropy.table import Table
from astropy.utils import iers
iers.conf.auto_download = False

from TTF import TTF_Query ##(DateObs, min_depth, Target)
##from TTF import TTF_Query_obj
from Utils import Get_Exp_Master, Check_AzEl, Get_Coo, Get_Times, Check_El

import warnings
warnings.simplefilter("ignore")

import textwrap
#################################################################
min_depth = 7.0
max_mag = 14.
Kou= EarthLocation(lat=57.036537*u.deg, \
                lon=59.545735*u.deg, height=290*u.m)

#################################################################
Local_Time = Time.now()
##Local_Time = Time('2022-08-02T12:28:01.0', format='fits')
##Test = TimeDelta(5*24*3600, format='sec')
##Local_Time = Local_Time + Test
TimeZone = TimeDelta(5*3600, format='sec')
Start = Time(np.round(Local_Time.jd), format='jd') - TimeZone

#################################################################
##get Sun and Moon
Plan = Table()
Step = TimeDelta(180, format='sec') 
Epochs = Time(np.arange(Start.jd, Start.jd+1., Step.jd), format='jd')
Plan['JD'] = Epochs

altazframe = AltAz(obstime=Plan['JD'], location=Kou)
sunaltaz = get_sun(Epochs).transform_to(altazframe)
Plan['Sun_El'] = sunaltaz.alt
Plan['Sun_Az'] = sunaltaz.az
moonaltaz = get_moon(Epochs).transform_to(altazframe)
Plan['Moon_El'] = moonaltaz.alt
Plan['Moon_Az'] = moonaltaz.az

Plan['Time'] = (Plan['JD'].jd - Plan['JD'][0].jd)*24. ##hours after midday

Zero = np.where(Plan['Sun_El']<=0)
Sunset = np.min(Plan['JD'][Zero].jd)
Sunrise = np.max(Plan['JD'][Zero].jd)
Title = 'Sunset: ' + Time(Sunset, format='jd').datetime.strftime('%Y-%m-%d %H:%M')
Title += ', '
Title += 'sunrise:' + Time(Sunrise, format='jd').datetime.strftime('%Y-%m-%d %H:%M')
Min_lim = np.min(Plan['Time'][Zero])
Max_lim = np.max(Plan['Time'][Zero])

##fig, axs = plt.subplots(2, 1, figsize=(8, 13), dpi=125)
fig = plt.figure(figsize=(8, 6), dpi=125)
# fig.add_axes([0, 0, 1, 1])
gs = gridspec.GridSpec(1, 1, height_ratios=[1])
axs0 = plt.subplot(gs[0])
# axs1 = plt.subplot(gs[1])

try:
    Night = np.where(Plan['Sun_El']<=-18)
    N_Start = np.min(Plan['Time'][Night])
    N_Stop = np.max(Plan['Time'][Night])
    axs[0].axvspan(N_Start, N_Stop, facecolor='k', alpha=0.3)
except:
    pass

try:
    Twilight = np.where(Plan['Sun_El']<-12)
    T_Start = np.min(Plan['Time'][Twilight])
    T_Stop = np.max(Plan['Time'][Twilight])
    axs0.axvspan(T_Start, T_Stop, facecolor='k', alpha=0.2)
    T_Start = Plan['JD'][Plan['Time']==T_Start]
    T_Stop = Plan['JD'][Plan['Time']==T_Stop]
except:
    pass

##axs0.set_position([0.125, 0.45, 0.8, 0.49])
axs0.hlines(-18, 0, 24, 'k', 'dashed', alpha=0.5)
axs0.hlines(-12, 0, 24, 'k', 'dashed', alpha=0.2)
axs0.hlines(30, 0, 24, 'r', 'dashed', alpha=0.3)

axs0.plot(Plan['Time'], Plan['Sun_El'], 'r-.', alpha=0.5, label='Sun')
axs0.plot(Plan['Time'], Plan['Moon_El'], 'b-.', alpha=0.5, label='Moon')

#################################################################
##axs1.set_position([0.125, 0.05, 0.8, 0.35])
collabel=('Target', 'Coord', 'Start(UTC)', 'Duration(h)', 'Depth(mmag)',  \
          'Vmag', 'Exp(s)', '2Moon(d)', 'Comments', 'Priority')
'''widths=([3/210., 5/210., 4/210., 2/210., 2/210., 1.5/210., 1.5/210., 2/210., 189/210.])
axs1.axis('tight')
axs1.axis('off')'''
Data = []
colors = []

##get TESS target list
DT = Local_Time.datetime.strftime('%m-%d-%Y')
Tess_List = TTF_Query(DT, min_depth, max_mag, '')
ascii.write(Tess_List, 'TTF_Info.txt', \
                fast_writer=False, overwrite=True, delimiter='\t', format='commented_header')
##Tess_List = TTF_Query_obj('207425167')
##print(Tess_List)
if len(Tess_List)>0:
    for row in range(len(Tess_List)):
        ##times
        t1 = Tess_List['jd_start'][row]+2450000 ##tr start
        t2 = Tess_List['jd_end'][row]+2450000   ##tr stop
        t0 = t1-0.75/24. ##obs start
        t3 = t2+0.75/24. ##obs stop
        ##check twilight
        t0, t1, t2, t3, Partial = Get_Times(t0, t1, t2, t3, T_Start, T_Stop)
        ##check partial observations and remove it from the list
        if Partial:
            continue
        ObsTime = Time(np.arange(t0, t3, Step.jd), format='jd')
        TrTime = Time(np.arange(t1, t2, Step.jd), format='jd')
        if (len(TrTime)>10) & (len(ObsTime)>10):
            #get position
            Ra, Dec = Get_Coo(Tess_List['coords(J2000)'][row])
            Eq = SkyCoord(ra = Ra*u.deg, dec = Dec*u.deg, frame='icrs')
            ObsAltAz = Eq.transform_to(AltAz(obstime = ObsTime, location = Kou))
            ObsAltAz = Table([ObsAltAz.obstime, ObsAltAz.alt, ObsAltAz.az], names=('obstime', 'El', 'Az'))
            TrAltAz = Eq.transform_to(AltAz(obstime = TrTime, location = Kou))
            TrAltAz = Table([TrAltAz.obstime, TrAltAz.alt, TrAltAz.az], names=('obstime', 'El', 'Az'))
            
            ##check elevation
            ObsAltAz, TrAltAz, Trim = Check_El(ObsAltAz, TrAltAz)
            ##check incimplite transits and remove from list
            if Trim > 0:
                continue
            ##checl moon
            Moon = get_moon(Time((t0+t3)/2., format='jd'))
            sep = Moon.separation(Eq).degree
            sep = np.round(sep, 2)
            
            colors.append(['w', 'w', 'w', 'w', 'w', 'w', 'w', 'w', 'w'])
##            if Check_AzEl(Tess_List['coords(J2000)'][row], Tess_List['az_start'][row]):
##                colors[-1][1] = 'lightpink' ##bad master angle
            Exp = Get_Exp_Master(Tess_List['V'][row])
            if Exp <=2:
                colors[-1][6] = 'lightpink' ##short exposure
            if Tess_List['V'][row] > 14:
                colors[-1][5] = 'lightpink' ##faint star
            if Tess_List['depth(mmag)'][row] < 6:
                colors[-1][4] = 'lightpink' ##min depth
            if sep<90:
                colors[-1][7] = 'lightpink' ##moon distance
            if Trim > 0:
                colors[-1][3] = 'lightpink' ##trimmed
                
            Begin = ObsAltAz['obstime'][0].datetime.strftime('%Y-%m-%d %H:%M:%S')
            Duration = round((ObsAltAz['obstime'][-1].jd - ObsAltAz['obstime'][0].jd)*24, 2)
            Data.append([Tess_List['Name'][row].replace('.', '_').replace(' ', ''),\
                         Tess_List['coords(J2000)'][row],\
                         Begin ,\
                         Duration,\
                         Tess_List['depth(mmag)'][row],\
                         Tess_List['V'][row],\
                         Get_Exp_Master(Tess_List['V'][row]),\
                         sep, \
                         Tess_List['comments'][row]
                        ])
            
            D = axs0.plot((ObsAltAz['obstime'].jd - Plan['JD'][0].jd)*24.,\
                        ObsAltAz['El'], \
                        linestyle = 'dotted', linewidth=2, label='_nolegend_')
            
            axs0.plot((TrAltAz['obstime'].jd - Plan['JD'][0].jd)*24.,\
                        TrAltAz['El'], linestyle = 'solid', \
                        label=Tess_List['Name'][row],\
                        color = D[-1].get_color(), linewidth=2)
##        colors[-1][0] = D[-1].get_color()
    '''
    the_table = axs1.table(cellText=Data, colLabels=collabel, \
                             colWidths=widths, cellColours = colors,\
                             loc='upper center')
    the_table.auto_set_font_size(False)
    the_table.set_fontsize(3)
    the_table.scale(2, 4) '''
    
    Title = Title + '<br>TESS transits with depth>'+str(min_depth)+'mmag and Vmag<'+str(max_mag)
else:
    Title = 'List is empty for some reasons. Maybe night starts/ends at nautical twilight.'

#################################################################
axs0.set_xticks(np.linspace(0, 23, 24))
axs0.set_xticklabels(['07', '08', '09', '10', '11', '12', '13', '14', '15', \
                                    '16', '17', '18', '19', '20', '21', '22', '23', '00', \
                                    '01', '02', '03', '04', '05', '06'])
axs0.tick_params(axis='both', labelsize=6, direction='in')
axs0.set_xlim(Min_lim, Max_lim)
axs0.set_xlabel('UTC, ' + Start.datetime.strftime('%Y-%m-%d'), fontsize=6)
axs0.set_ylim(-20, None)
axs0.set_ylabel('Elevation (deg)', fontsize=6)
axs0.legend(loc='lower center', fontsize=5, bbox_to_anchor=(0.5, 1), ncol=6)
axs0.grid()


##plt.show()

# fig.suptitle(Title, fontsize=8)

##Name = 'PlanT.pdf'
##if os.path.isfile(Name):
##   os.remove(Name)
##plt.savefig(Name)
##shutil.move(Name, '/var/www/htdocs/PlanT.pdf')

plt.savefig("PlanPlot.svg", format="svg", transparent=True, bbox_inches='tight')
# plt.show()


def rateObservations(observationInfo = []):
    rating = 0
    rating = +observationInfo[collabel.index("2Moon(d)")] * moonLuminatedPercent/100 + observationInfo[collabel.index("Depth(mmag)")] + observationInfo[collabel.index("Exp(s)")] - observationInfo[collabel.index("Vmag")]
    # print(observationInfo[collabel.index("Depth(mmag)")], observationInfo[collabel.index("Vmag")], observationInfo[collabel.index("Exp(s)")], observationInfo[collabel.index("2Moon(d)")], rating)
    # observationInfo.append(rating)

    return rating


'''
collabelsForDisplay = {"Target": "Name", "Coord": "coords(J2000)", "Start(UTC)": ""}

def getHTMLTableFromAstropy(astropyTable):
    columnHeaders = astropyTable.colnames

    print(columnHeaders, collabel)
    table = '<table>'
    table += '<tr>'
    for header in columnHeaders:
        if header in collabelsForDisplay.values():
            table += '<th>' + str(header) + '</th>'
    table += '</tr>'

    for row in range(len(astropyTable)):
        table += '<tr>'
        for col in columnHeaders:
            if col in collabelsForDisplay.values():
                table += '<th>' + str(astropyTable[row][col]) + '</th>'
        table += '</tr>'
    
    table += '</table>'
    return table
'''

def getHTMLTable(columnHeaders = [''], rows = [['']]):
    colOfTarget = 0
    colOfComments = 0

    table = '<table>'
    table += '<tr>'
    for i in range(len(columnHeaders)):
        if columnHeaders[i] == "Name":
            colOfTarget = i
        if columnHeaders[i] == "Comments":
            colOfComments = i
            continue
        table += '<th>' + str(columnHeaders[i]) + '</th>'
    table += '</tr>'

    for row in rows:
        table += '<tr>'
        for i in range(len(row)):
            if i == colOfComments:
                continue
            if i == colOfTarget:
                table += '<th><a class="targetName" href="comments.html?com='+ urllib.parse.quote_plus(str(row[colOfComments])) +'&targ='+ urllib.parse.quote_plus(str(row[i])) +'">' + str(row[i]) + '</a></th>'
            else:
                table += '<th>' +str(row[i]) +'</th>'
        table += '</tr>'
    
    table += '</table>'
    return table 


# print(getHTMLTable(collabel, Data))
# print(repr(Tess_List.keys()))
# Tess_List.show_in_browser(jsviewer=True)



date = ephem.now()
nnm = ephem.next_new_moon(date)
pnm = ephem.previous_new_moon(date)
moonPhase = (date-pnm)/(nnm-pnm)

moon = ephem.Moon()
moon.compute()
moonLuminatedPercent = moon.phase

Data.sort(key=rateObservations, reverse=True)
# collabel.append("Priority")
for i in range(len(Data)):
    Data[i].append(str(i+1))




pageTemplateFile = open("pageTemplate.html", "r")
pageTemplate = pageTemplateFile.read()
renderedPage = open("index.html", "w")
# print(pageTemplate)

renderedPage.write(Template(pageTemplate).substitute(Title = Title, mainTable = getHTMLTable(collabel, Data), MoonPhase = moonPhase, MoonLuminatedPercent = int(moonLuminatedPercent)))

pageTemplateFile.close()
renderedPage.close()