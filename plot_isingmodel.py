import pandas as pd
import sys
import matplotlib.pyplot as plt

lx = int(sys.argv[1])
model = str(sys.argv[2])

df = pd.read_csv(f"{model}_Size{lx}.csv")
temp_data = df['Temperature'].tolist()
mag_data = df['Magnetisation'].tolist()
sus_data = df['Susceptibility'].tolist()
eng_data = df['Energy'].tolist()
shc_data = df['HeatCapacity'].tolist()
mag_err = df['mag_error'].tolist()
sus_err = df['sus_error'].tolist()
eng_err = df['eng_error'].tolist()
shc_err = df['shc_error'].tolist()

fig, axs = plt.subplots(2, 2)
if model == "Glauber": fig.suptitle('Glauber, Size 50x50')
else: fig.suptitle('Kawasaki, Size 50x50')

axs[0, 0].errorbar(temp_data, mag_data, marker = "o", color = "dodgerblue", linestyle = '', markersize = 4, yerr = mag_err)
#axs[0, 0].set_title('Magnetisation')
axs[0, 0].set_ylabel('Average Magnetisation')

axs[0, 1].errorbar(temp_data, eng_data, marker = "o", color = "dodgerblue", linestyle = '', markersize = 4, yerr = eng_err)
#axs[0, 1].set_title('Energy')
axs[0, 1].set_ylabel('Average Energy')

axs[1, 0].errorbar(temp_data, sus_data, marker = "o", color = "dodgerblue", linestyle = '', markersize = 4, yerr = sus_err)
#axs[1, 0].set_title('Susceptibility')
axs[1, 0].set_ylabel('Susceptibility')

axs[1, 1].errorbar(temp_data, shc_data, marker = "o", color = "dodgerblue", linestyle = '', markersize = 4, yerr = shc_err)
#axs[1, 1].set_title('Specific Heat Capacity')
axs[1, 1].set_ylabel('Specific Heat Capacity')

for ax in axs.flat:
    ax.set(xlabel='Temperature')

fig.tight_layout(pad=1.0)

plt.show()