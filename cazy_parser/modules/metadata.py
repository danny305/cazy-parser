from pathlib import Path
import logging

import pandas as pd


log = logging.getLogger("cazylog")

def dump_metadata(metadata: list[dict], output_f: Path):
    """
    Save the metadata to a file.

    Parameters
    ----------
    metadata : list
        List of metadata dictionaries.
    output_f : str
        Path to the output file.
    """
    df = pd.DataFrame(metadata)
    df.to_csv(output_f, index=False)

    log.info(f"Saved cazy metadata: {output_f}")
