import data_handling
import matplotlib.pyplot as plt
dict_20x20 	= data_handling.open_saved(20, 4000)
beta_list 	= dict_20x20['beta']
mag_list 	= dict_20x20['M']
sus_list 	= dict_20x20['chi']

fig, ax = plt.subplots(2,figsize = (10,10))
for axis in ax:
	axis.set_xlabel(r"$\beta / J$", fontsize = 20)
ax[0].set_ylabel("Magnetization (a.u.)", fontsize = 20, labelpad = 20)
ax[1].set_ylabel("Susceptibility (a.u.)", fontsize = 20, labelpad = 20)
ax[0].plot(beta_list, mag_list, 'xb', markersize=6)
ax[1].plot(beta_list, sus_list, 'xr', markersize=6)
ax[0].axvline(x=.41, ymin=0, ymax=140,linewidth =.4,c='k',linestyle='--')
ax[1].axvline(x=.41, ymin=0, ymax=140,linewidth =.4,c='k',linestyle='--')
fig.savefig("Results/20x20_4000.png", dpi = 300)
