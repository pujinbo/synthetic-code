###################################################################
# Replication file for Cattaneo, Feng, Palomba, and Titiunik (2022)
###################################################################

########################################
# Load SCPI_PKG package
import pandas
import numpy
import random
import os
from warnings import filterwarnings
from plotnine import ggtitle, ggsave, theme, element_blank
from copy import deepcopy

from scpi_pkg.scdata import scdata
from scpi_pkg.scdataMulti import scdataMulti
from scpi_pkg.scest import scest
from scpi_pkg.scpi import scpi
from scpi_pkg.scplot import scplot
from scpi_pkg.scplotMulti import scplotMulti

os.chdir('YOUR_PATH_HERE')
filterwarnings('ignore')
sims = 200
cores = 1

########################################
# Load database
data = pandas.read_stata('final_data.dta')
data = data.loc[data['restricted2'] == 1]
data['lgdp'] = numpy.log(data['rgdppp'])
data['year'] = pandas.DatetimeIndex(data['year']).year  # from datetime64 to int

########################################
# Analyze separately each continent

for method in ["L1-L2"]:

    for cont in ["Asia", "Africa", "Europe", "North America", "South America"]:

        units = data.loc[(data['continent'] == cont) & (data['treated'] == 1) &
                         (data['trDate'] <= 1992)]['countryname'].unique().tolist()

        if cont == "Europe":
            units = [un for un in units if un != "Slovenia"]

        ###############################################
        # unit-time treatment effect
        ###############################################
        print(method + ' - ' + cont + ': unit-time')
        aux = scdataMulti(df=data,
                          id_var='countryname',
                          treatment_var='liberalization',
                          outcome_var='lgdp',
                          time_var='year',
                          features={'treated': ['lgdp', 'lsc']},
                          constant=True,
                          cointegrated_data=True,
                          cov_adj={'treated': ['constant', 'trend']},
                          post_est=10,
                          units_est=units,
                          effect='unit-time',
                          anticipation=1)

        res = scpi(aux, w_constr={'name': 'simplex'}, sims=sims,
                   u_order=1, u_lags=1, cores=cores)
        p = scplotMulti(res, ptype='series', scales='free_y', ncols=2)

        p = (p + theme(axis_text_y=element_blank(),
                       axis_ticks_major_y=element_blank()))
        plotname = cont.replace(' ', '_') + '_unittime_' + method + '.png'
        ggsave(filename=plotname, plot=p)

        ###############################################
        # average unit post-treatment effect
        ###############################################
        print(method + ' - ' + cont + ': unit')
        aux = scdataMulti(df=data,
                          id_var='countryname',
                          treatment_var='liberalization',
                          outcome_var='lgdp',
                          time_var='year',
                          features={'treated': ['lgdp', 'lsc']},
                          constant=True,
                          cointegrated_data=True,
                          cov_adj={'treated': ['constant', 'trend']},
                          post_est=10,
                          units_est=units,
                          effect='unit')

        res = scpi(aux, w_constr={'name': 'simplex'}, sims=sims,
                   u_order=1, u_lags=1, cores=cores)
        p = scplotMulti(res, ptype='series')
        plotname = cont.replace(' ', '_') + '_unit_' + method + '.png'
        ggsave(filename=plotname, plot=p)

        ###############################################
        # average post-treatment effect on the treated
        ###############################################
        print(method + ' - ' + cont + ': time')
        aux = scdataMulti(df=data,
                          id_var='countryname',
                          treatment_var='liberalization',
                          outcome_var='lgdp',
                          time_var='year',
                          features={'treated': ['lgdp', 'lsc']},
                          constant=True,
                          cointegrated_data=True,
                          cov_adj={'treated': ['constant', 'trend']},
                          post_est=10,
                          units_est=units,
                          effect='time')

        res = scpi(aux, w_constr={'name': 'simplex'}, sims=sims,
                   u_order=1, u_lags=1, cores=cores)

        p = scplotMulti(res, ptype='series')
        plotname = cont.replace(' ', '_') + '_time_' + method + '.png'
        ggsave(filename=plotname, plot=p)
