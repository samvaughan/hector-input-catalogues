"""
Make a redshift catalogue from all of our Hector Redshift Survey observations. Only keep spectra which have a redshift probability greater than 0.95.
The output catalogue is saved as a parquet file
"""

import pandas as pd
import utils
from pathlib import Path
from tqdm import tqdm


if __name__ == "__main__":
    smk = snakemake  # noqa
    all_files = list(smk.input.redshift_observations)

    print(
        f"Loading {len(all_files)} redshift observation files and keeping observations with PROB > 0.95..."
    )
    df = pd.DataFrame()
    for filename in tqdm(all_files):
        filename = Path(filename)
        tmp_df = utils.load_FITS_table_in_pandas(filename)

        field_name = filename.stem.split("_")[0]
        assert (
            field_name in smk.params.all_field_names
        ), f"We seem to have a field which we didn't expect: {field_name}"
        # Make sure the field is the first column
        tmp_df.insert(0, "FIELD", field_name)
        df = pd.concat((df, tmp_df))

    df = df.loc[df.PROB > 0.95]
    df.to_parquet(smk.output.hrs_redshift_catalogue)
    print(f"...Done! {len(df)} redshifts pass the quality cut")
