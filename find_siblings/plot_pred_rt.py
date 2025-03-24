import os
import pandas as pd
from koinapy import Koina

path = "CSV_mouse_rat"

# Read all files from path in a single dataframe
# Ignore empty files
files = os.listdir(path)
dfs = []
for file in files:
    if os.path.getsize(f"{path}/{file}"):
        df = pd.read_csv(f"{path}/{file}", header=None)
        dfs.append(df)

# Concatenate all dataframes
df = pd.concat(dfs, ignore_index=True)

# Change column names
df.columns = ["Peptide_1", "Peptide_2"]

# Create input files for Koina
inputs_1 = pd.DataFrame()
inputs_2 = pd.DataFrame()
inputs_1["peptide_sequences"] = df["Peptide_1"]
inputs_2["peptide_sequences"] = df["Peptide_2"]

# Predict the retention time using Koina
# Can also use different models such as Prosit_2019_irt, AlphaPeptDeep_rt_generic, ...
model = Koina("Deeplc_hela_hf", "koina.wilhelmlab.org:443")
predictions_1 = model.predict(inputs_1)
predictions_2 = model.predict(inputs_2)

# Add the predictions to the dataframe
df["Peptide_1_RT"] = predictions_1["irt"]
df["Peptide_2_RT"] = predictions_2["irt"]
df["length"] = df["Peptide_1"].str.len()

# Calculate the normalized delta RT
df["delta_RT"] = df["Peptide_2_RT"] - df["Peptide_1_RT"]
df["absolute_delta_RT"] = df["delta_RT"].abs()
df["average_RT"] = (df["Peptide_1_RT"] + df["Peptide_2_RT"]) / 2
df["delta_RT_normalized"] = df["absolute_delta_RT"] / df["average_RT"]

max = df["absolute_delta_RT"].max()
df[df["delta_RT_normalized"] > 0.3][
    ["Peptide_1", "Peptide_2", "Peptide_1_RT", "Peptide_2_RT"]
]

# Plot length on x-axis and delta_RT on y-axis
# Add the mean delta RT for each length as a red dot
import matplotlib.pyplot as plt

plt.scatter(df["length"], df["delta_RT_normalized"])
plt.scatter(
    df.groupby("length")["delta_RT_normalized"].mean().index,
    df.groupby("length")["delta_RT_normalized"].mean().values,
    color="red",
)
plt.xlabel("Peptide length")
plt.ylabel("Delta RT/Average RT")
plt.title("Delta RT vs peptide length")
plt.savefig("delta_rt_vs_length.png")
plt.close()
