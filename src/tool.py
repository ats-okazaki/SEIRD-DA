import numpy as np
import pandas as pd
import param as par

def getdata(country,key) :

    key2 = {'lon':'Longitude (average)','lat':'Latitude (average)'}
    #
    # read file
    #
    df1 = pd.read_csv('/home/okazaki/tmp/covid19/data/org/owid-covid-data.csv', index_col=0)
    df2 = pd.read_csv('/home/okazaki/tmp/covid19/data/org/countries_codes_and_coordinates.csv', index_col=2)

    #
    # extract data
    #
    if key == 'population' or key == 'location' :
        dat = df1.at[country,key][0]
    elif key == 'lon' or key == 'lat' :
        try :
            dat = df2.at[' \"'+country+'\"',key2[key]]
            if type(dat) is not str :
                dat = dat[0]
            dat = float(dat.replace('"',''))
        except KeyError :
            print('Not found ',country,' --skipped')
            dat = -999.
    return (dat)
#------------------------------------------------------------------------
def setobs_hdx(country,date) :
    #
    # initialize
    #
    '''
    yo[0] : newly confirmed case
    yo[1] : newly recovered case
    yo[2] : newly death case
    '''
    yo = np.empty(3,dtype=float)
    year = int(date.astype(object).year) - 2000
    mon  = int(date.astype(object).month)
    day  = int(date.astype(object).day)
    today = str(mon)+ '/' +str(day)+ '/' +str(year)

    #
    # read file
    #
    df_i = pd.read_csv('/home/okazaki/tmp/covid19/data/org/time_series_covid19_confirmed_global_iso3_regions.csv')
    df_r = pd.read_csv('/home/okazaki/tmp/covid19/data/org/time_series_covid19_recovered_global_iso3_regions.csv')
    df_d = pd.read_csv('/home/okazaki/tmp/covid19/data/org/time_series_covid19_deaths_global_iso3_regions.csv')

    #
    # extract data
    #
    try :
        tot1 = df_i[ df_i['ISO 3166-1 Alpha 3-Codes'] == country ][today]
        yo[0] = np.sum(tot1)
    except KeyError :
        yo[0] = par.undef

    try :
        tot1 = df_r[ df_r['ISO 3166-1 Alpha 3-Codes'] == country ][today]
        yo[1] = np.sum(tot1) 
    except KeyError :
        yo[1] = par.undef

    try :
        tot1 = df_d[ df_d['ISO 3166-1 Alpha 3-Codes'] == country ][today]
        yo[2] = np.sum(tot1) 
    except KeyError :
        yo[2] = par.undef

    return ( yo )
#------------------------------------------------------------------------
def setobs_oxf(country,date) :

    yo = np.empty(2,dtype=float)
    #
    # read file
    #
    df = pd.read_csv('/home/okazaki/tmp/covid19/data/org/owid-covid-data.csv', index_col=2)

    #
    # extract data
    #
    try :
        yo[0] = df[ df['iso_code'] == country ].at[ str(date),'new_cases' ]
        yo[1] = df[ df['iso_code'] == country ].at[ str(date),'new_deaths' ]
    except KeyError :
        yo[:] = par.undef
    return (yo)
    
