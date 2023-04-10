from astropy.time import Time
import requests
from astropy.io import ascii
import re

def create_TTF_entry_url(date_obs, min_depth, max_mag, target):  
    url = 'https://astro.swarthmore.edu/telescope/tess-secure/' + \
          'print_eclipses.cgi?observatory_string=57.036537%3B' + \
          '59.545735%3BAsia%2FYekaterinburg%3BKourovka+' + \
          'Observatory+0.4m+%28x2%29%3BKourovka+Obs+0.4m+' + \
          '%28x2%29&use_utc=1&observatory_latitude=57.036537' + \
          '&observatory_longitude=59.545735&timezone=UTC' + \
          '&start_date={}&days_to_print=1&days_in_past=0' + \
          '&minimum_start_elevation=30' + \
          '&and_vs_or=or&minimum_end_elevation=30' + \
          '&minimum_ha=-12&maximum_ha=12&baseline_hrs=0' + \
          '&show_unc=1&maximum_priority=3&minimum_depth={}' + \
          '&maximum_V_mag={}&target_string={}' + \
          '&lco_only=0&single_object=0&ra=&dec=&epoch=' + \
          '&period=&duration=&target=&show_ephemeris=0' + \
          '&print_html=2&twilight=-12&max_airmass=2.4'
          
    return(url.format(date_obs, min_depth, max_mag, target))

def TTF_Query(DateObs, min_depth, max_mag, Target):
    TTF_URL = create_TTF_entry_url(DateObs, min_depth, max_mag, Target)
    ssn = requests.session()
    ssn.auth = ('tess_nda_observer', 'F1nd_TE$S_PlaNets!')
    req = ssn.get(TTF_URL)
    if req.status_code != 200:
        message = None
    else:
##        print(req.text)
##        Cleared = req.text
##        lines = req.text.splitlines()
##        for line in lines:
##            print(line)
##            print()
##            Index = [i for i, ltr in enumerate(line) if ltr == '"']
##            if len(Index)>0:
##                Comment = '"'+line[Index[0]+1:Index[-1]].replace('"', 'arcsec')+'"'
##                line = line[:Index[0]]+Comment
##            Cleared =  Cleared + line + '\n'
##        Clear = re.sub(r'[^\x00-\x7f]',r'', req.text)
        message = ascii.read(re.sub(r'[^\x00-\x7f]',r'', req.text), delimiter = ',', format = 'commented_header', quotechar = '"')
##        message = ascii.read(req.text, delimiter = ',', format = 'commented_header', quotechar = '"') 
##        re.sub(r'[^\x00-\x7f]',r'', message)
        
    print(message[0][25])
    return message


def TTF_Query_obj(DateObs, Target):
    url = 'https://astro.swarthmore.edu/telescope/tess-secure/' + \
          'print_eclipses.cgi?observatory_string=57.036537%3B' + \
          '59.545735%3BAsia%2FYekaterinburg%3BKourovka+' + \
          'Observatory+0.4m+%28x2%29%3BKourovka+Obs+0.4m+' + \
          '%28x2%29&use_utc=1&observatory_latitude=57.036537' + \
          '&observatory_longitude=59.545735&timezone=UTC' + \
          '&start_date={}&days_to_print=1&days_in_past=0' + \
          '&minimum_start_elevation=30' + \
          '&and_vs_or=or&minimum_end_elevation=30' + \
          '&minimum_ha=-12&maximum_ha=12&baseline_hrs=0' + \
          '&show_unc=1' + \
          '&target_string={}' + \
          '&lco_only=0&single_object=0&ra=&dec=&epoch=' + \
          '&period=&duration=&target=&show_ephemeris=0' + \
          '&print_html=2&twilight=-12&max_airmass=2.4'
    TTF_URL = url.format(DateObs, Target)
    ssn = requests.session()
    ssn.auth = ('tess_nda_observer', 'F1nd_TE$S_PlaNets!')
    req = ssn.get(TTF_URL)
    if req.status_code != 200:
        message = None
    else:
##        print(req.text)
##        Cleared = req.text
##        lines = req.text.splitlines()
##        for line in lines:
##            print(line)
##            print()
##            Index = [i for i, ltr in enumerate(line) if ltr == '"']
##            if len(Index)>0:
##                Comment = '"'+line[Index[0]+1:Index[-1]].replace('"', 'arcsec')+'"'
##                line = line[:Index[0]]+Comment
##            Cleared =  Cleared + line + '\n'
        message = ascii.read(req.text, delimiter = ',', format = 'commented_header', quotechar = '"') #, \
                             #delimiter = ',', format = 'commented_header')
    # print(message)
    return message
##DateObs = '2022-08-02T12:28:01.0'
##T = Time(DateObs, format='fits', out_subfmt='longdate')
##DT = T.datetime.strftime('%m-%d-%Y')
##Target = 'TIC207425167_01'
##Target = Target[0:3]+'+'+Target[3:].replace('_', '.')
##print(Target)
####
##print(TTF_Query_obj(DT, Target))
##print(Q.info())
##

