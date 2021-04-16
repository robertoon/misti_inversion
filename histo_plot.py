import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import seaborn
from math import sqrt
fontname = 'Helvetica'
sb_dark = seaborn.dark_palette('skyblue', 8, reverse=True)
#seaborn.set(palette=sb_dark)
np.seterr(divide='ignore', invalid='ignore')

fig, ax = plt.subplots(nrows=4, ncols=2, figsize=(7.48, 9.05))
#bins = 20

df_post = pd.read_csv('filtered_posterior_1_5phi.csv')
df_prior = pd.read_csv('prior.csv')

########### Column height ###################################
post = df_post["column_height"]
prior = df_prior["column_height"]
ch_po, ch_pr = post, prior

n, min_max, mean, var, skew, kurt = stats.describe(post)
std = sqrt(var)
confidence = stats.t.interval(0.95,len(post)-1,loc=mean, scale=std/sqrt(len(post)))

height, bins, patches = ax[0,0].hist(ch_pr, bins = 80, density=True, lw=0.1, ec='k', color='#668A05', label="prior")
height, bins, patches = ax[0,0].hist(ch_po, bins = 80, density=True, lw=0.1, ec='k', color='#056588', label="posterior")

ax[0,0].axvline(x= mean - (2 * std), ymin=0, ymax=5, color='black',linestyle='--',linewidth=0.8)
ax[0,0].axvline(x= mean + (2 * std), ymin=0, ymax=5, color='black',linestyle='--',linewidth=0.8)
ax[0,0].yaxis.set_visible(False)
#ax[0,0].text(6500, 0.00017, 'a')
ax[0,0].xaxis.set_tick_params(labelsize=8)
ax[0,0].xaxis.get_major_formatter().set_powerlimits((0, 1))
ax[0,0].xaxis.offsetText.set_fontsize(8)
ax[0,0].set_title('Umbrella height (m.a.s.l)', fontname=fontname, fontsize=8)

########### Disk radius ###################################
post = df_post["ellipse_major_axis"]
prior = df_prior["ellipse_major_axis"]
ema_po, ema_pr = post, prior

n, min_max, mean, var, skew, kurt = stats.describe(post)
std = sqrt(var)
confidence = stats.t.interval(0.95,len(post)-1,loc=mean, scale=std/sqrt(len(post)))

height, bins, patches = ax[0,1].hist(ema_pr, bins = 60, density=True, lw=0.1, ec='k', color='#CEE4ED', label="prior")
height, bins, patches = ax[0,1].hist(ema_po, bins = 60, density=True, lw=0.1, ec='k', color='#34768E', label="posterior")

ax[0,1].axvline(x= mean - (2 * std), ymin=0, ymax=5, color='black',linestyle='--',linewidth=0.8)
ax[0,1].axvline(x= mean + (2 * std), ymin=0, ymax=5, color='black',linestyle='--',linewidth=0.8)
ax[0,1].yaxis.set_visible(False)
ax[0,1].xaxis.get_major_formatter().set_powerlimits((0, 1))
ax[0,1].xaxis.set_tick_params(labelsize=8)
ax[0,1].xaxis.offsetText.set_fontsize(8)
ax[0,1].set_title('Umbrella major axis (m)', fontname=fontname, fontsize=8)

post = df_post["ellipse_major_axis"]
prior = df_prior["ellipse_major_axis"]
emi_po, emi_pr = post, prior

n, min_max, mean, var, skew, kurt = stats.describe(post)
std = sqrt(var)
confidence = stats.t.interval(0.95,len(post)-1,loc=mean, scale=std/sqrt(len(post)))

height, bins, patches = ax[1,0].hist(emi_pr, bins = 60, density=True, lw=0.1, ec='k', color='#CEE4ED', label="prior")
height, bins, patches = ax[1,0].hist(emi_po, bins = 60, density=True, lw=0.1, ec='k', color='#34768E', label="posterior")

ax[1,0].axvline(x= mean - (2 * std), ymin=0, ymax=5, color='black',linestyle='--',linewidth=0.8)
ax[1,0].axvline(x= mean + (2 * std), ymin=0, ymax=5, color='black',linestyle='--',linewidth=0.8)
ax[1,0].yaxis.set_visible(False)
ax[1,0].xaxis.get_major_formatter().set_powerlimits((0, 1))
ax[1,0].xaxis.set_tick_params(labelsize=8)
ax[1,0].xaxis.offsetText.set_fontsize(8)
ax[1,0].set_title('Umbrella minor axis (m)', fontname=fontname, fontsize=8)

# ########### Erupted mass ###################################
post = df_post["total_erupted_mass"]
prior = df_prior["total_erupted_mass"]
em_po, em_pr = post, prior

n, min_max, mean, var, skew, kurt = stats.describe(post)
std = sqrt(var)
confidence = stats.t.interval(0.95,len(post)-1,loc=mean, scale=std/sqrt(len(post)))

height, bins, patches = ax[1,1].hist(em_pr, bins = 60, density=True, lw=0.1, ec='k', color='#CEE4ED', label="prior")
height, bins, patches = ax[1,1].hist(em_po, bins = 60, density=True, lw=0.1, ec='k', color='#34768E', label="posterior")

ax[1,1].axvline(x= mean - (2 * std), ymin=0, ymax=5, color='black',linestyle='--',linewidth=0.8)
ax[1,1].axvline(x= mean + (2 * std), ymin=0, ymax=5, color='black',linestyle='--',linewidth=0.8)
ax[1,1].yaxis.set_visible(False)
ax[1,1].xaxis.get_major_formatter().set_powerlimits((0, 1))
ax[1,1].xaxis.set_tick_params(labelsize=8)
ax[1,1].xaxis.offsetText.set_fontsize(8)
ax[1,1].set_title('Erupted mass (kg)', fontname=fontname, fontsize=8)

# ########### TGSD mean ###################################
post = df_post["tgsd_mean"]
prior = df_prior["tgsd_mean"]
gs_po, gs_pr = post, prior

n, min_max, mean, var, skew, kurt = stats.describe(post)
std = sqrt(var)
confidence = stats.t.interval(0.95,len(post)-1,loc=mean, scale=std/sqrt(len(post)))

height, bins, patches = ax[2,0].hist(gs_pr, bins = 60, density=True, lw=0.1, ec='k', color='#CEE4ED', label="prior")
height, bins, patches = ax[2,0].hist(gs_po, bins = 60, density=True, lw=0.1, ec='k', color='#34768E', label="posterior")

ax[2,0].axvline(x= mean - (2 * std), ymin=0, ymax=5, color='black',linestyle='--',linewidth=0.8)
ax[2,0].axvline(x= mean + (2 * std), ymin=0, ymax=5, color='black',linestyle='--',linewidth=0.8)
ax[2,0].yaxis.set_visible(False)
ax[2,0].xaxis.set_tick_params(labelsize=8)
ax[2,0].set_title('TGSD mean ($\phi$)', fontname=fontname, fontsize=8)

# ########### Diffusion ###################################
post = df_post["diffusion_coef"]
prior = df_prior["diffusion_coef"]
df_po, df_pr = post, prior

n, min_max, mean, var, skew, kurt = stats.describe(post)
std = sqrt(var)
confidence = stats.t.interval(0.95,len(post)-1,loc=mean, scale=std/sqrt(len(post)))

height, bins, patches = ax[2,1].hist(df_pr, bins = 60, density=True, lw=0.1, ec='k', color='#CEE4ED', label="prior")
height, bins, patches = ax[2,1].hist(df_po, bins = 60, density=True, lw=0.1, ec='k', color='#34768E', label="posterior")

ax[2,1].axvline(x= mean - (2 * std), ymin=0, ymax=5, color='black',linestyle='--',linewidth=0.8)
ax[2,1].axvline(x= mean + (2 * std), ymin=0, ymax=5, color='black',linestyle='--',linewidth=0.8)
ax[2,1].yaxis.set_visible(False)
ax[2,1].xaxis.get_major_formatter().set_powerlimits((0, 1))
ax[2,1].xaxis.set_tick_params(labelsize=8)
ax[2,1].xaxis.offsetText.set_fontsize(8)
ax[2,1].set_title('Diffusion ($m^{2} s^{-1}$)', fontname=fontname, fontsize=8)

# ########### Wind speed ###################################
post = df_post["wind_speed"]
prior = df_prior["wind_speed"]
ws_po, ws_pr = post, prior

n, min_max, mean, var, skew, kurt = stats.describe(post)
std = sqrt(var)
confidence = stats.t.interval(0.95,len(post)-1,loc=mean, scale=std/sqrt(len(post)))

height, bins, patches = ax[3,0].hist(ws_pr, bins = 60, density=True, lw=0.1, ec='k', color='#CEE4ED', label="prior")
height, bins, patches = ax[3,0].hist(ws_po, bins = 60, density=True, lw=0.1, ec='k', color='#34768E', label="posterior")

ax[3,0].axvline(x= mean - (2 * std), ymin=0, ymax=5, color='black',linestyle='--',linewidth=0.8)
ax[3,0].axvline(x= mean + (2 * std), ymin=0, ymax=5, color='black',linestyle='--',linewidth=0.8)
ax[3,0].yaxis.set_visible(False)
ax[3,0].xaxis.set_tick_params(labelsize=8)
ax[3,0].set_title('Wind speed ($m s^{-1}$)', fontname=fontname, fontsize=8)

# ########### Wind direction ###################################
post = df_post["wind_direction"]
prior = df_prior["wind_direction"]
wd_po, wd_pr = post, prior

n, min_max, mean, var, skew, kurt = stats.describe(post)
std = sqrt(var)
confidence = stats.t.interval(0.95,len(post)-1,loc=mean, scale=std/sqrt(len(post)))

height, bins, patches = ax[3,1].hist(wd_pr, bins = 60, density=True, lw=0.1, ec='k', color='#CEE4ED', label="prior")
height, bins, patches = ax[3,1].hist(wd_po, bins = 60, density=True, lw=0.1, ec='k', color='#34768E', label="posterior")

ax[3,1].axvline(x= mean - (2 * std), ymin=0, ymax=5, color='black',linestyle='--',linewidth=0.8)
ax[3,1].axvline(x= mean + (2 * std), ymin=0, ymax=5, color='black',linestyle='--',linewidth=0.8)
ax[3,1].yaxis.set_visible(False)
ax[3,1].xaxis.set_tick_params(labelsize=8)
ax[3,1].set_title('Wind direction ($^{0}$)', fontname=fontname, fontsize=8)
###################################################################
#fig.suptitle('Best_phi * 1.1')
fig.tight_layout(pad=0.5) 
plt.show() 
#plt.savefig('Best_phi * 1,5.png', dpi=1800, facecolor='w', edgecolor='k')