import pandas as pd

smk = snakemake  # noqa
waves_S = pd.read_parquet(smk.input.waves_S_cat)
waves_N = pd.read_parquet(smk.input.waves_N_cat)
cluster_cat = pd.read_parquet(smk.input.cluster_cat)

waves_S["Field"] = "WAVES_S"
waves_N["Field"] = "WAVES_N"
cluster_cat["Field"] = "HectorClusters"

combined_catalogue = pd.concat((waves_N, waves_S, cluster_cat))
combined_catalogue.to_csv(smk.output.master_cat, index=False)
