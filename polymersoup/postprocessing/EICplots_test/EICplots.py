import pandas as pd
from pandas.io.json import json_normalize
import json
import matplotlib.pyplot as plt
import seaborn as sns

# filepath
EICs_json = 'ECK2-11_HGSL5_1in10_PRE_FILTERED_MS1_EICs.json'

# read EICs.json
with open(EICs_json, 'r') as EICs_json_data:
    EICs_dict = json.load(EICs_json_data)

# create pandas dataframe from json data
EICs_normalized_data = pd.DataFrame.from_dict(json_normalize(EICs_dict), orient='columns')

# change dataframe layout (sequence becomes indexed variable)
EICs_data = pd.melt(EICs_normalized_data)

# count number of rows
nrows = len(EICs_data)

# create new columns for retention times and intensities
EICs_data["retention times"] = EICs_data["value"]
EICs_data["intensities"] = EICs_data["value"]

# separate retention times and intensities to different rows
for i in range(nrows):
    rt = list(zip(*EICs_data["value"][i]))[0]
    EICs_data["retention times"][i] = rt
    intensity = list(zip(*EICs_data["value"][i]))[1]
    EICs_data["intensities"][i] = intensity

# remove 'value' column, rename 'variable' to sequence
EICs_data = EICs_data.drop(columns = 'value')
EICs_data.rename(columns = {'variable':'sequence'}, inplace=True)

# create EIC plot for each sequence (retention time vs intensity)
sns.set_style("ticks")
sns.set_palette("Set2")
for i in range(nrows):
    plt.close()
    plt.xlabel("Retention Time")
    plt.ylabel("Intensity")
    plt.title(EICs_data["sequence"][i])
    sns.lineplot((EICs_data["retention times"][i]),  (EICs_data["intensities"][i]), data=EICs_data)
    # plt.savefig("%s_%s.png" % (i , EICs_data["sequence"][i]))