#!/usr/bin/env python
# coding: utf-8

# In[1]:


#1 Dependencies and Setup
import matplotlib.pyplot as plt
import pandas as pd
import numpy as mp
import scipy.stats as st


# In[2]:


#1 Study data files
mouse_metadata_path = "data/Mouse_metadata.csv"
study_results_path = "data/Study_results.csv"


# In[3]:


#1 Read the mouse data and the study results
mouse_metadata = pd.read_csv(mouse_metadata_path)
study_results = pd.read_csv(study_results_path)


# In[4]:


#1 
mouse_metadata


# In[5]:


#1
study_results


# In[6]:


#1 Combine the data into a single DataFrame
all_data = pd.merge(study_results, mouse_metadata, how ='left', on = ['Mouse ID'])


# In[7]:


#1 Display the data table for preview
all_data.head()


# In[8]:


#2 Checking the number of mice.
all_data['Mouse ID'].nunique()


# In[9]:


#3 Our data should be uniquely identified by Mouse ID and Timepoint
# Get the duplicate mice by ID number that shows up for Mouse ID and Timepoint. 
duplicated = all_data[all_data.duplicated(subset=['Mouse ID', 'Timepoint'], keep= False)]['Mouse ID'].unique()
duplicated


# In[10]:


#4 Optional: Get all the data for the duplicate mouse ID. 
all_data.loc[all_data['Mouse ID']=='g989']


# In[11]:


#5 Create a clean DataFrame by dropping the duplicate mouse by its ID.
all_data = all_data[all_data['Mouse ID'] != 'g989']


# In[12]:


#5 Create a clean DataFrame by dropping the duplicate mouse by its ID.
all_data = all_data[all_data['Mouse ID'] != 'g989']
all_data['Mouse ID'].nunique()
all_data.head()


# In[13]:


#6 Checking the number of mice in the clean DataFrame.
all_data['Mouse ID'].nunique()


# ## Summary Statistics

# In[14]:


#7 Generate a summary statistics table of mean, median, variance, standard deviation, and SEM of the tumor volume for each regimen

# Use groupby and summary statistical methods to calculate the following properties of each drug regimen: 
# mean, median, variance, standard deviation, and SEM of the tumor volume. 
# Assemble the resulting series into a single summary DataFrame.


# In[15]:


#7 mean
mean= all_data.groupby('Drug Regimen')['Tumor Volume (mm3)'].mean()
mean


# In[16]:


#7 median
median= all_data.groupby('Drug Regimen')['Tumor Volume (mm3)'].median()
median


# In[17]:


#7 variance
variance= all_data.groupby('Drug Regimen')['Tumor Volume (mm3)'].var()
variance


# In[18]:


#7 standard deviation
std= all_data.groupby('Drug Regimen')['Tumor Volume (mm3)'].std()
std


# In[19]:


#7 SEM
sem= all_data.groupby('Drug Regimen')['Tumor Volume (mm3)'].sem()


# In[20]:


#7 Assemble the resulting series into a single summary DataFrame.

statistics = pd.DataFrame({
    'Mean Tumor': mean,
    'Median Tumor Size': median,
    'Variance': variance,
    "Standard Deviation": std,
    "SEM": sem
})
statistics


# In[21]:


#8 A more advanced method to generate a summary statistics table of mean, median, variance, standard deviation,
# and SEM of the tumor volume for each regimen (only one method is required in the solution)
# Using the aggregation method, produce the same summary statistics in a single line

all_data.groupby('Drug Regimen').agg({'Tumor Volume (mm3)': ['mean', 'median', 'var', 'std', 'sem']})


# In[22]:


#9 Generate a bar plot showing the total number of rows (Mouse ID/Timepoints) for each drug regimen using Pandas.
drug_counts = all_data.groupby('Drug Regimen')['Timepoint'].count()
drug_counts = drug_counts.sort_values(ascending=False)

ax = drug_counts.plot(kind="bar", figsize=(7, 5), color="blue", fontsize=10)
ax.set_xlabel("Drug Regimen", fontsize=12)
ax.set_ylabel("# of Observed Mouse Timepoints", fontsize=12)
plt.show()


# In[23]:


#10 Generate a bar plot showing the total number of rows (Mouse ID/Timepoints) for each drug regimen using pyplot.
drug_counts = all_data.groupby('Drug Regimen')['Timepoint'].count()
drug_counts = drug_counts.reset_index().rename(columns={'Timepoint': '# of Observed Mouse Timepoints'})
drug_counts = drug_counts.sort_values(by='# of Observed Mouse Timepoints', ascending=False)

drug_counts


# In[24]:


#10 Generate a bar plot showing the total number of rows (Mouse ID/Timepoints) for each drug regimen using pyplot.
x_value = drug_counts['Drug Regimen']
label = drug_counts['# of Observed Mouse Timepoints']

plt.bar(x_value, label)

plt.xlabel('Drug Regimen')
plt.xticks(rotation='vertical')
plt.ylabel('# of Observed Mouse Timepoints')

plt.show()


# In[25]:


#11 Generate a pie plot showing the distribution of female versus male mice using Pandas
sex_counts = all_data['Sex'].value_counts()

sex_counts.plot(kind='pie', autopct='%1.1f%%')


# In[26]:


#12 Generate a pie plot showing the distribution of female versus male mice using pyplot
sexes = all_data['Sex'].value_counts()
sexes


# In[27]:


#12 Generate a pie plot showing the distribution of female versus male mice using pyplot
labels = ['Males', 'Females']
sizes = [958, 922]
colors = ['blue', 'orange']

plt.pie(sizes, labels = labels, autopct = '%1.1f%%')
plt.show


# ## Quartiles, Outliers and Boxplots

# In[28]:


#13 Calculate the final tumor volume of each mouse across four of the treatment regimens:  
#13 Capomulin, Ramicane, Infubinol, and Ceftamin

#13 Start by getting the last (greatest) timepoint for each mouse


#13 Merge this group df with the original DataFrame to get the tumor volume at the last timepoint


# In[29]:


#13 Start by getting the last (greatest) timepoint for each mouse

last_time = all_data.groupby(['Drug Regimen','Mouse ID'])['Timepoint'].max()
last_time = last_time.reset_index().rename(columns={'Timepoint': 'Last Timepoint'})

last_time


# In[30]:


#13 Merge this group df with the original DataFrame to get the tumor volume at the last timepoint
last_time= pd.merge(last_time, all_data, how= 'left')
last_time


# In[31]:


#13 Capomulin, Ramicane, Infubinol, and Ceftamin
treatments = ['Capomulin', 'Ramicane', 'Infubinol', 'Ceftamin']
drugs_filtered = last_time[last_time['Drug Regimen'].isin(treatments)]
drugs_filtered


# In[32]:


#14 Put treatments into a list for for loop (and later for plot labels)
#14  Create empty list to fill with tumor vol data (for plotting)
treatments = ['Capomulin', 'Ramicane', 'Infubinol', 'Ceftamin']
tumor_volume = []

for treatment in treatments:
    last_call = drugs_filtered[drugs_filtered['Drug Regimen']== treatment]
        
    time_point = last_call.groupby('Mouse ID')['Timepoint'].max()

    final_tumor_volume = pd.merge(time_point, last_call, on=['Mouse ID', 'Timepoint'])['Tumor Volume (mm3)']

    
    tumor_volume.append(final_tumor_volume)

final_tumor_volume = pd.DataFrame(tumor_volume).transpose()
final_tumor_volume.columns = treatments

final_tumor_volume


# In[33]:


#14 Put treatments into a list for for loop (and later for plot labels)
treatments = ['Capomulin', 'Ramicane', 'Infubinol', 'Ceftamin']
outliers = []

for drug in treatments:
    # tumor volume data for the current drug
    volumes = final_tumor_volume[drug]
    
    # Get quartiles and IQR
    q1 = volumes.quantile(0.25)
    q3 = volumes.quantile(0.75)
    iqr = q3 - q1
    
    # Get lower and upper bounds for outliers
    lower_bound = q1 - (1.5 * iqr)
    upper_bound = q3 + (1.5 * iqr)
    
    # Find  potential outliers
    outliers_mask = (volumes < lower_bound) | (volumes > upper_bound)
    drug_outliers = volumes[outliers_mask]
    
    # Append results with the list of outliers
    outliers.append({
        'Treatment': drug,
        'First Quartile': q1,
        'Third Quartile': q3,
        'IQR': iqr,
        'Lower Bound': lower_bound,
        'Upper Bound': upper_bound,
        'Minimum': volumes.min(),
        'Maximum': volumes.max(),
        'Median': volumes.median(),
        'Outliers': drug_outliers
    })

outliers = pd.DataFrame(outliers)
outliers


# In[34]:


#14
outliers_grouped = outliers.groupby('Treatment')['Outliers'].agg(list)
outliers_grouped


# In[35]:


#15 Generate a box plot that shows the distrubution of the tumor volume for each treatment group.
cap_data = final_tumor_volume['Capomulin']
ram_data = final_tumor_volume['Ramicane']
inf_data = final_tumor_volume['Infubinol']
cef_data = final_tumor_volume['Ceftamin']

data = [cap_data, ram_data, inf_data, cef_data]

fig, ax = plt.subplots()
ax.boxplot(data)

ax.set_xticklabels(['Capomulin', 'Ramicane', 'Infubinol', 'Ceftamin'])

ax.set_ylabel('Final Tumor Volume (mm3)')

ax.set_title('Boxplots of Final Tumor Volume by Treatment')

flierprops = dict(marker='D', markerfacecolor='b', markersize=10)

plt.boxplot(data, flierprops=flierprops)
plt.show()


# In[36]:


#16 Generate a line plot of tumor volume vs. time point for a single mouse treated with Capomulin
single_mouse = all_data.loc[(all_data['Drug Regimen'] == 'Capomulin') & (all_data['Mouse ID'] == 'l509')]

plt.plot(single_mouse['Timepoint'], single_mouse['Tumor Volume (mm3)'])

plt.xlabel('Timepoints (Days)')
plt.ylabel('Tumor Volume (mm3)')
plt.title('Capomulin Treatment Treatment of Mouse l509')
plt.show()


# In[37]:


#17 Generate a scatter plot of mouse weight vs. the average observed tumor volume for the entire Capomulin regimen
capomulin_data = all_data.loc[all_data['Drug Regimen'] == 'Capomulin']
cap_avg = capomulin_data.groupby(['Mouse ID']).mean()

plt.scatter(cap_avg['Weight (g)'],cap_avg['Tumor Volume (mm3)'])

plt.xlabel('Weight (g)')
plt.ylabel('Average Tumor Volume (mm3)')

plt.show()


# In[38]:


#18 Calculate the correlation coefficientfor mouse weight and average observed tumor volume for the entire Capomulin regimen
corr=cap_avg['Weight (g)'].corr(cap_avg['Tumor Volume (mm3)'])
corr= round(corr, 2)

print(f"The correlation between mouse weight and average tumor volume is {corr}")


# In[39]:


#18 linear regression model for mouse weight and average observed tumor volume for the entire Capomulin regimen
import scipy.stats as st

lineregress=st.linregress(cap_avg['Weight (g)'], cap_avg['Tumor Volume (mm3)'])
lineregress


# In[40]:


#18
x_values = cap_avg['Weight (g)']
y_values = cap_avg['Tumor Volume (mm3)']

slope, intercept, rvalue, pvalue, stderr = st.linregress(cap_avg['Weight (g)'], cap_avg['Tumor Volume (mm3)'])


# In[41]:


#18 Plot the linear regression model on top of the previous scatter plot.

y_predict = slope *(cap_avg['Weight (g)'])+intercept

plt.scatter(cap_avg['Weight (g)'],cap_avg['Tumor Volume (mm3)'])

plt.plot(cap_avg['Weight (g)'], y_predict, 'r'.format(slope,intercept))

plt.xlabel('Weight(g)')
plt.ylabel('Average Tumore Volume (mm3)')

plt.show()


# In[42]:


#14 an example of how much time a For Loop can save you... accidently did this first haha

#Ramicane
iqr_cap = final_tumor_volume['Capomulin'].describe()
iqr_cap

q1_cap = iqr_cap['25%']
median_cap = iqr_cap['50%']
q3_cap = iqr_cap['75%']

iqr_cap = q3_cap - q1_cap

lower_bound_cap = q1_cap - (iqr_cap*1.5)
upper_bound_cap = q3_cap + (iqr_cap*1.5)

min_cap = final_tumor_volume[treatments[0]].min()
max_cap = final_tumor_volume[treatments[0]].max()

#Ramicane
iqr_ram = final_tumor_volume['Ramicane'].describe()
iqr_ram

q1_ram = iqr_ram['25%']
median_ram = iqr_ram['50%']
q3_ram = iqr_ram['75%']

iqr_ram = q3_ram - q1_ram

lower_bound_ram = q1_ram - (iqr_ram*1.5)
upper_bound_ram = q3_cap + (iqr_ram*1.5)

min_ram = final_tumor_volume[treatments[1]].min()
max_ram = final_tumor_volume[treatments[1]].max()

#Infubinol
iqr_inf = final_tumor_volume['Infubinol'].describe()
iqr_inf

q1_inf = iqr_inf['25%']
median_inf = iqr_inf['50%']
q3_inf = iqr_inf['75%']

iqr_inf = q3_inf - q1_inf

lower_bound_inf = q1_inf - (iqr_inf*1.5)
upper_bound_inf = q3_inf + (iqr_inf*1.5)

min_inf = final_tumor_volume[treatments[2]].min()
max_inf = final_tumor_volume[treatments[2]].max()

#Ceftamin
iqr_cef = final_tumor_volume['Ceftamin'].describe()
iqr_cef

q1_cef = iqr_cef['25%']
median_cef = iqr_cef['50%']
q3_cef = iqr_cef['75%']

iqr_cef = q3_cef - q1_cef

lower_bound_cef = q1_cef - (iqr_cef*1.5)
upper_bound_cef = q3_cef + (iqr_cef*1.5)

min_cef = final_tumor_volume[treatments[3]].min()
max_cef = final_tumor_volume[treatments[3]].max()

outliers = pd.DataFrame ({'Treatment': treatments,
                         'First Quartile': [q1_cap, q1_ram, q1_inf, q1_cef],
                         'Third Quartile': [q3_cap, q3_ram, q3_inf, q3_cef],
                         'IQR': [iqr_cap, iqr_ram, iqr_inf, iqr_cef],
                         'Lower Bound': [lower_bound_cap, lower_bound_ram, lower_bound_inf, lower_bound_cef],
                         'Upper Bound': [upper_bound_cap, upper_bound_ram, upper_bound_inf, upper_bound_cef],
                         'Minimum': [min_cap, min_ram, min_inf, min_cef],
                         'Maximum': [max_cap, max_ram, max_inf, max_cef],
                         'Median': [median_cap, median_ram, median_inf, median_cef]
                        })

outliers

