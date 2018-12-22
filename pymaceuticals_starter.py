#!/usr/bin/env python
# coding: utf-8

# In[1]:


# Dependencies and Setup
get_ipython().run_line_magic('matplotlib', 'inline')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.stats import sem

# Hide warning messages in notebook
import warnings
warnings.filterwarnings('ignore')

# File to Load 
mouse_drug_data_to_load = "data/mouse_drug_data.csv"
clinical_trial_data_to_load = "data/clinicaltrial_data.csv"

# Read the Mouse and Drug Data and the Clinical Trial Data
mouse_drug_df = pd.read_csv("data/mouse_drug_data.csv")
clinical_trial_df = pd.read_csv("data/clinicaltrial_data.csv")

# Combine the data into a single dataset
pyma_combine_df = pd.merge(clinical_trial_df,mouse_drug_df, how ='left' , on = 'Mouse ID')
pyma_combine_df=pyma_combine_df.drop_duplicates()

# Display the data table for preview
pyma_combine_df.head()


# ## Tumor Response to Treatment

# In[2]:


# Store the Mean Tumor Volume Data Grouped by Drug and Timepoint 
tumor_volume_df = pyma_combine_df.groupby(['Drug','Timepoint']).mean()['Tumor Volume (mm3)']
tumor_volume_df=tumor_volume_df.reset_index()

# Convert to DataFrame
tumor_volume_df = pd.DataFrame(tumor_volume_df)

# Preview DataFrame
tumor_volume_df


# In[3]:


# Store the Standard Error of Tumor Volumes Grouped by Drug and Timepoint
tumor_volume_sem = pyma_combine_df.groupby(['Drug','Timepoint']).sem()['Tumor Volume (mm3)']
tumor_volume_sem = tumor_volume_sem.reset_index()

# Convert to DataFrame
tumor_volume_sem = pd.DataFrame(tumor_volume_sem)

# Preview DataFrame
tumor_volume_sem.head()


# In[4]:


# Minor Data Munging to Re-Format the Data Frames
tumor_volume_df = pd.pivot_table(tumor_volume_df, index='Timepoint', columns='Drug', values='Tumor Volume (mm3)', aggfunc = np.mean)

# Preview that Reformatting worked
tumor_volume_df.head()


# In[5]:


# Generate the Plot (with Error Bars)
Timepoint = tumor_volume_df.index
plt.errorbar(Timepoint, tumor_volume_df["Capomulin"], yerr=tumor_volume_df["Capomulin"].sem(), color="r", marker="o", markersize=5, linestyle="--", linewidth=0.50)
plt.errorbar(Timepoint, tumor_volume_df["Infubinol"], yerr=tumor_volume_df["Infubinol"].sem(), color="b", marker="^", markersize=5, linestyle="--", linewidth=0.50)
plt.errorbar(Timepoint, tumor_volume_df["Ketapril"], yerr=tumor_volume_df["Ketapril"].sem(), color="g", marker="s", markersize=5, linestyle="--", linewidth=0.50)
plt.errorbar(Timepoint, tumor_volume_df["Placebo"], yerr=tumor_volume_df["Placebo"].sem(), color="k", marker="d", markersize=5, linestyle="--", linewidth=0.50)

plt.title("Tumor Response to Treatment")
plt.xlabel("Time (Days)")
plt.ylabel("Tumor Volume (mm3)")
lim = (0,max(Timepoint))
plt.legend(loc = 'best', frameon=True)
plt.grid(True)

# Save the Figure
plt.savefig("Analysis/Fig1.png")

# Show the Figure
plt.show()


# ## Metastatic Response to Treatment

# In[6]:


# Store the Mean Met. Site Data Grouped by Drug and Timepoint 
met_site_mean = pyma_combine_df.groupby(["Drug","Timepoint"]).mean()["Metastatic Sites"]

# Convert to DataFrame
met_site_mean = pd.DataFrame(met_site_mean)

# Preview DataFrame
met_site_mean.head()


# In[7]:


# Store the Standard Error associated with Met. Sites Grouped by Drug and Timepoint 
met_site_std = pyma_combine_df.groupby(["Drug","Timepoint"]).sem()["Metastatic Sites"]

# Convert to DataFrame
met_site_std = pd.DataFrame(met_site_std)

# Preview DataFrame
met_site_std.head()


# In[8]:


# Minor Data Munging to Re-Format the Data Frames
met_site_mean = met_site_mean.reset_index()
met_site_mean_pivot = met_site_mean.pivot(index = "Timepoint",columns ="Drug")["Metastatic Sites"]

met_site_std = met_site_std.reset_index()
met_site_std_pivot = met_site_std.pivot(index = "Timepoint",columns ="Drug")["Metastatic Sites"]

# Preview that Reformatting worked
met_site_mean_pivot.head()


# In[9]:


# Generate the Plot (with Error Bars)
Metastatic = met_site_mean_pivot.index
plt.errorbar(Metastatic, met_site_mean_pivot["Capomulin"], yerr= met_site_std_pivot["Capomulin"], color="r", marker="o", markersize=5, linestyle="--", linewidth=0.5)
plt.errorbar(Metastatic, met_site_mean_pivot["Infubinol"], yerr= met_site_std_pivot["Infubinol"], color="b", marker="^", markersize=5, linestyle="--", linewidth=0.5)
plt.errorbar(Metastatic, met_site_mean_pivot["Ketapril"], yerr= met_site_std_pivot["Ketapril"], color="g", marker="s", markersize=5, linestyle="--", linewidth=0.5)
plt.errorbar(Metastatic, met_site_mean_pivot["Placebo"], yerr= met_site_std_pivot["Placebo"], color="k", marker="d", markersize=5, linestyle="--", linewidth=0.5)

plt.title("Metastatic Spread During Treatment")
plt.xlabel("Treatment Duration (Days)")
plt.ylabel("Met. Sites")
lim = (0,max(Timepoint))
plt.legend(loc = 'upper left', frameon=True)
plt.grid(True)

# Save the Figure
plt.savefig("Analysis/Fig2.png")

# Show the Figure
plt.show()


# ## Survival Rates

# In[10]:


# Store the Count of Mice Grouped by Drug and Timepoint (W can pass any metric)
survival_rate= pyma_combine_df.groupby(["Drug","Timepoint"]).count()["Tumor Volume (mm3)"]

# Convert to DataFrame
survival_rate = pd.DataFrame({"Mouse Count" : survival_rate})
survival_rate =survival_rate.reset_index()

# Preview DataFrame
survival_rate.head()


# In[11]:


# Minor Data Munging to Re-Format the Data Frames
survival_rate = survival_rate.reset_index()
survival_rate_pivot = survival_rate.pivot(index = "Timepoint",columns = "Drug")["Mouse Count"] 

# Preview the Data Frame
survival_rate_pivot.head()


# In[12]:


# Generate the Plot (Accounting for percentages)
plt.plot(100 * survival_rate_pivot["Capomulin"] /25,'ro',linestyle ="--", markersize =5, linewidth =0.50)
plt.plot(100 * survival_rate_pivot["Infubinol"] /25,'b^',linestyle ="--", markersize =5, linewidth =0.50)
plt.plot(100 * survival_rate_pivot["Ketapril"] /25,'gs',linestyle ="--", markersize =5, linewidth =0.50)
plt.plot(100 * survival_rate_pivot["Placebo"] /25,'kd',linestyle ="--", markersize =5, linewidth =0.50)

plt.title("Survival During Treatment")
plt.xlabel("Time (Days)")
plt.ylabel("Survival Rate (%)")
plt.legend(loc = "bottom left", fontsize ="small")
plt.grid(True)

# Save the Figure
plt.savefig("Analysis/Fig3.png")

# Show the Figure
plt.show()


# ## Summary Bar Graph

# In[13]:


# Calculate the percent changes for each drug
tumor_pct_change =100*(tumor_volume_df.iloc[-1] - tumor_volume_df.iloc[0]) / tumor_volume_df.iloc[0]

# Display the data to confirm
tumor_pct_change


# In[14]:


# Store all Relevant Percent Changes into a Tuple
perct_changes = (tumor_pct_change["Capomulin"],
               tumor_pct_change["Infubinol"],
               tumor_pct_change["Ketapril"],
               tumor_pct_change["Placebo"])

# Splice the data between passing and failing drugs
fig, ax = plt.subplots()
ind = np.arange(len(perct_changes))
width = 1
rectsPass = ax.bar(ind[0], perct_changes[0], width, color='green')
rectsFail = ax.bar(ind[1:], perct_changes[1:], width, color='red')

# Orient widths. Add labels, tick marks, etc.
plt.ylabel('% Tumor Volume Change')
plt.title('Tumor Change Over 45 Day Treatment')
plt.ylim([-30,70])
plt.grid(True)
plt.xticks(ind + 0.5)

ax.set_xticklabels(('Capomulin', 'Infubinol', 'Ketapril', 'Placebo'))


# Use functions to label the percentages of changes
def autolabelFail(rects):
    for rect in rects:
        height = rect.get_height()
        ax.text(rect.get_x() + rect.get_width()/2., 3,
                '%d%%' % int(height),
                ha='center', va='bottom', color="white")

def autolabelPass(rects):
    for rect in rects:
        height = rect.get_height()
        ax.text(rect.get_x() + rect.get_width()/2., -8,
                '-%d%% ' % int(height),
                ha='center', va='bottom', color="white")

# Call functions to implement the function calls
autolabelPass(rectsPass)
autolabelFail(rectsFail)

# Save the Figure
plt.savefig("Analysis/Fig4.png")

# Show the Figure
fig.show()


# In[ ]:




