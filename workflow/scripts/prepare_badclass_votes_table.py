"""
Using the community votes from the hector target selection app (https://hts.datacentral.org.au/hector), 
prepare a table of the _probability_ that each target falls into a given 'badclass'
"""
import json
import numpy as np
import pandas as pd
import scipy.special as sp
from tqdm import tqdm
from collections import defaultdict

smk = snakemake  # noqa


def get_votes_per_classification(votes):
    """Get the number of votes each class has in a JSON file

    Args:
        votes: Votes from a JSON file

    Returns:
        array-like: The counts for each of the 9 different categories.
    """
    n = np.zeros(9)
    for i in range(9):
        n[i] = np.count_nonzero(votes == i)

    return n


with open(smk.input.json_votes, "r") as f:
    data = json.load(f)

N_galaxies = len(data)  # 2832

category_probabilities = np.zeros((N_galaxies - 1, 9))
galaxy_catids = [d["target"] for d in data[1:]]
N_classifiers = np.zeros(N_galaxies - 1)

# Count classifers by user
user_classifications = defaultdict(list)
for d in data[1:]:
    classifications = d["badclass"]
    classifications.pop("default", None)
    for key, value in classifications.items():
        value["target"] = d["target"]
        user_classifications[key].append(value)

# Find the number of classifications per person
for key, value in user_classifications.items():
    print(f"{key}: {len(value)}")

# A dictionary for saving each result
data_dict = {}
for i in tqdm(range(1, N_galaxies)):
    galaxy = pd.DataFrame(data[i]["badclass"]).T
    # Remove the "default" class
    galaxy = galaxy.loc[galaxy.index != "default"]

    # Find the number of classifiers
    N_classifiers[i - 1] = len(galaxy)
    votes = pd.to_numeric(galaxy.value).values
    alpha = get_votes_per_classification(votes) + 0.5

    # Get the probabilities that each category is the true consensus vote.
    # Using a helpful document prepared for Caro!
    # see file:///Users/samvaughan/Desktop/Mark_Donoghoe_report.pdf,
    log_p_k = np.zeros(9)
    for k in range(0, 9):
        for j in range(9):
            if j == k:
                pass
            else:
                log_p_k[k] += np.log(sp.betainc(alpha[j], alpha[k], 0.5))
    p_k = np.exp(log_p_k)
    category_probabilities[i - 1, :] = p_k

    # Make a dictionary from these values
    galaxy_data_dict = {}
    galaxy_data_dict["ID"] = data[i]["target"]
    galaxy_data_dict["p_badclass_0"] = category_probabilities[i - 1, 0]
    galaxy_data_dict["p_badclass_1"] = category_probabilities[i - 1, 1]
    galaxy_data_dict["p_badclass_2"] = category_probabilities[i - 1, 2]
    galaxy_data_dict["p_badclass_3"] = category_probabilities[i - 1, 3]
    galaxy_data_dict["p_badclass_4"] = category_probabilities[i - 1, 4]
    galaxy_data_dict["p_badclass_5"] = category_probabilities[i - 1, 5]
    galaxy_data_dict["p_badclass_6"] = category_probabilities[i - 1, 6]
    galaxy_data_dict["p_badclass_7"] = category_probabilities[i - 1, 7]
    galaxy_data_dict["p_badclass_8"] = category_probabilities[i - 1, 8]

    data_dict[data[i]["target"]] = galaxy_data_dict

master_df = pd.DataFrame(data_dict).T

# Add in the most probable bad class, along with its probability
master_df["BAD_CLASS"] = master_df.apply(
    lambda x: np.argmax(x.filter(regex="p_badclass_*")), axis=1
)
master_df["p_BAD_CLASS"] = master_df.apply(
    lambda x: np.max(x.filter(regex="p_badclass_*")), axis=1
)
master_df.to_csv(smk.output.badclass_probabilities, index=False)
