import numpy as np

import pandas as pd

from sklearn.decomposition import PCA

import plotly.express as px

import matplotlib.pyplot as plt

from matplotlib.legend_handler import HandlerLine2D

# Define the LegendColorHandler class
class LegendColorHandler(HandlerLine2D):
    def create_artists(self, legend, orig_handle,
                       xdescent, ydescent, width, height, fontsize, trans):
        line = plt.Line2D([xdescent, xdescent + width], [ydescent, ydescent], color=orig_handle)
        return [line]

def change_results(data):
    FT = data['FT']
    PHASE_MARGIN = data['PHASE_MARGIN']
    results = data['all_satured']
    for i in range(len(results)):
        if results[i] == 1 and FT[i] >= 13.7e06 and PHASE_MARGIN[i] >= 60:
            results[i] = 2
    
    data['all_satured'] = results
    
            
def color (results):
    c = []
    for i in range(len(results)):
        if results[i] == 0:
            c.append('red')
        elif results[i] == 1:
            c.append('blue')
        elif results[i] == 2:
            c.append('green')
    return c
    
def plot_pca(data, n) :
    pca = PCA(n_components=n)
    principalComponents = pca.fit_transform(data)
    pc1 = principalComponents[:,0]
    pc2 = principalComponents[:,1]
    c = color(data['all_satured'])
    
    plt.scatter(pc1, pc2, c=c)
    
    # Add a legend
    #legend_labels = ['Fail in DC', 'Fail in AC', 'Succed in AC and DC']  # update the labels to match your data
    #legend_colors = ['red', 'blue', 'green']  # update the colors to match your plot
    #plt.legend(legend_labels, title='Circuit performance', bbox_to_anchor=(0., -0.3, 1., .102), loc=3, ncol=3, mode="expand", borderaxespad=0., fancybox=True, shadow=True, title_fontsize=12, facecolor='#FFFFFF', edgecolor='#333333', handler_map={legend_labels[0]: LegendColorHandler(legend_colors[0]), legend_labels[1]: LegendColorHandler(legend_colors[1]), legend_labels[2]: LegendColorHandler(legend_colors[2])})
 
    plt.show()
    


datas = np.genfromtxt('results/output_file_optimisation.txt', delimiter=' ', skip_header=1)

headers = ['gmid1', 'gmid3', 'gmid5', 'gmid6', 'gmid7', 'gmid8', 'all_satured', 'AV0', 'FT', 'PHASE_MARGIN']

dataset = pd.DataFrame(datas, columns=headers)

change_results(dataset)
print(dataset)

fig = plot_pca(dataset, 2)

filter_succes = dataset[dataset.all_satured == 2]
filter_DC_fail = dataset[dataset.all_satured == 0]
filter_AC_fail = dataset[dataset.all_satured == 1]

scale_y = 1 
scale_x = 1e-6
plt.plot(filter_succes['FT']*scale_x, filter_succes['PHASE_MARGIN']*scale_y, 'o', c='green', label='Succed in AC and DC')
plt.plot(filter_DC_fail['FT']*scale_x, filter_DC_fail['PHASE_MARGIN']*scale_y, 'o', c='red', label='Fail in DC')
plt.plot(filter_AC_fail['FT']*scale_x, filter_AC_fail['PHASE_MARGIN']*scale_y, 'o', c='blue', label='Fail in AC')
plt.xlabel('FT [MHz]')
plt.ylabel('PHASE MARGIN [°]')
plt.legend()
plt.show()

scale_y = 1 
scale_x = 1e-6
# Plot the first data series using the left y-axis
plt.plot(filter_succes['AV0']*scale_y, filter_succes['FT']*scale_x, 'o', c='purple', label='FT [MHz]')

plt.xlabel('Gain [dB]')
plt.ylabel('FT [MHz]')  # left y-axis label
plt.legend()

# Create a new y-axis on the right side of the plot
ax2 = plt.twinx()

# Plot the second data series using the right y-axis
ax2.plot(filter_succes['AV0']*scale_y, filter_succes['PHASE_MARGIN']*scale_y, '*', c='teal', label='Phase Margin [°]')

# Add labels to the y-axes
ax2.set_ylabel('AV0 [dB]')  # right y-axis label


# Add a legend
plt.legend()

# Show the plot
plt.show()


import seaborn as sns

# Compute the correlation matrix
corr = dataset.corr()

# Create a heatmap visualization
sns.heatmap(corr, cmap="PuOr", annot=True)

# Show the plot
plt.show()





