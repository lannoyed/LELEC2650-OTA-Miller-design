import numpy as np

import pandas as pd

import matplotlib.pyplot as plt

datas = np.genfromtxt('result.txt', delimiter=' ', skip_header=1)

headers = ['gmid16','gmid17', 'gmid12','B12-13', 'Flag_Correct_Saturation', 'Flag_Correct_Design', 'Flag_result_is_physical', 'Flag_mathematic_integrity']



dataset = pd.DataFrame(datas, columns=headers)





# Create a figure with 3 subplots
fig, axs = plt.subplots(1, 4, figsize=(16, 4))

# Count the number of occurrences of each unique value in the 'Flag_Correct_Saturation' column
counts = dataset['Flag_Correct_Saturation'].value_counts()

# Create a bar plot of the counts in the first subplot
counts.plot.bar(ax=axs[0])

# Count the number of occurrences of each unique value in the 'Flag_Correct_Design' column
counts = dataset['Flag_Correct_Design'].value_counts()

# Create a bar plot of the counts in the second subplot
counts.plot.bar(ax=axs[1])

# Count the number of occurrences of each unique value in the 'Flag_result_is_physical' column
counts = dataset['Flag_result_is_physical'].value_counts()

# Create a bar plot of the counts in the third subplot
counts.plot.bar(ax=axs[2])

counts = dataset['Flag_mathematic_integrity'].value_counts()

counts.plot.bar(ax=axs[3])

axs[0].set_title('Flag_Correct_Saturation')
axs[1].set_title('Flag_Correct_Design')
axs[2].set_title('Flag_result_is_physical')
axs[3].set_title('Flag_mathematic_integrity')


# Show the plot
plt.show()


import seaborn as sns

# Compute the correlation matrix
corr = dataset.corr()

# Create a heatmap visualization
sns.heatmap(corr, cmap="PuOr", annot=True)

# Show the plot
plt.show()

dataset.boxplot(column=['B12-13'], by='Flag_result_is_physical', figsize=(10, 6))
plt.show()