# -*- coding: utf-8 -*-
"""
Created on Fri Jan  5 10:52:52 2024

@author: Carlo Vignali

@def: This is primarily a notes page with information copied from ChatGPT


In Matplotlib, you have a wide range of elements in a plot that you can adjust. 
Below is a list of common elements you can customize in a Matplotlib plot, along 
with sample code snippets for how to make these adjustments.

"""

import matplotlib.pyplot as plt

### 1. Title of the Plot

plt.title("Plot Title", fontsize=14, color='blue')


### 2. Axis Labels

plt.xlabel("X-axis Label", fontsize=12)
plt.ylabel("Y-axis Label", fontsize=12)


### 3. Line Style and Thickness in Line Plots

plt.plot(x, y, color='green', linestyle='--', linewidth=2)


### 4. Axis Limits

plt.xlim(0, 10)
plt.ylim(-5, 5)


### 5. Tick Marks

plt.xticks([0, 2, 4, 6, 8, 10], rotation=45)
plt.yticks([-5, 0, 5], rotation=45)


### 6. Grid Lines

plt.grid(True, which='both', linestyle='-.', linewidth=0.5)


### 7. Legends

plt.legend(["Line 1"], loc='upper right')


### 8. Removing or Customizing Spines (Borders)

ax = plt.gca()
ax.spines['top'].set_visible(False)   # Remove top spine
ax.spines['right'].set_color('red')   # Change color of right spine
ax.spines['left'].set_linewidth(2)    # Change width of left spine


### 9. Figure Size

plt.figure(figsize=(8, 6))


### 10. Subplots

fig, ax = plt.subplots(2, 2)  # 2x2 grid of subplots


### 11. Scatter Plot Customization

plt.scatter(x, y, marker='o', color='purple', s=50)  # 's' is size of markers


### 12. Histogram Customization

plt.hist(data, bins=20, color='gray', alpha=0.7)


### 13. Boxplot Customization

plt.boxplot(data, patch_artist=True, notch=True)


### 14. Error Bars

plt.errorbar(x, y, yerr=error_array, fmt='o', color='black', ecolor='lightgray', elinewidth=3, capsize=0)


### 15. Text and Annotations

plt.text(5, 0, "Sample Text", fontsize=12, color='green')
plt.annotate('Annotation', xy=(3, 4), xytext=(5, 5), arrowprops=dict(facecolor='black'))


### 16. Customizing Tick Parameters

plt.tick_params(axis='both', which='major', labelsize=10, direction='inout', length=6)


### 17. Logarithmic Scale

plt.xscale('log')
plt.yscale('log')


### 18. Saving Plots

plt.savefig('plot.png', dpi=300, bbox_inches='tight')


