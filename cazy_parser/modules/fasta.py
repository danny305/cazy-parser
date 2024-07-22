from pathlib import Path
from time import time
import logging
import os

from Bio import Entrez, SeqIO
from tqdm import tqdm

Entrez.email = "your-email-here@example.org"

log = logging.getLogger("cazylog")


def download_fastas(id_list: list[str]) -> list[str]:
    """
    Download fasta files from NCBI.

    Parameters
    ----------
    id_list : list
        List of genbank ids.

    Returns
    -------
    fasta_list : list
        List of fasta strings.

    """
    start_t = time()
    log.info(f"Dowloading {len(id_list)} fasta sequences...")
    fasta_list = []
    handle = Entrez.efetch(db="protein", id=id_list, rettype="fasta", retmode="text")
    log.info(f"Entrez api call: {handle.url}")
    for seq_record in tqdm(SeqIO.parse(handle, "fasta"), total=len(id_list)):
        fasta_str = (
            f">{seq_record.description}{os.linesep}{seq_record.seq}{os.linesep*2}"
        )
        fasta_list.append(fasta_str)

    log.info(f"Downloaded {len(fasta_list)} fasta sequences in {time()-start_t:.2f} seconds")

    if len(fasta_list) != len(id_list):
        log.warning(
            f"Only {len(fasta_list)} fasta sequences were downloaded out of {len(id_list)}" 
        )
        
    return fasta_list


def dump_fastas(id_list: list[str], output_f: Path):
    """
    Save the fasta strings to a file.

    Parameters
    ----------
    id_list : list
        List of genbank ids.
    output_f : str
        Path to the output file.

    """
    assert isinstance(output_f, Path), f"{output_f} is of type {type(output_f)}"
    fasta_list = download_fastas(id_list)
    log.info(f"Dumping fasta sequences to file {output_f.resolve()}")
    with output_f.open("w") as fh:
        for fasta in fasta_list:
            fh.write(fasta)


def dump_id_list(id_list: list[str], output_f: Path) -> None:
    """Save the id list to a file."""
    assert isinstance(output_f, Path), f"{output_f} is of type {type(output_f)}"
    log.info(f"Dumping genbank ids to file {output_f.resolve()}")
    with output_f.open("w") as fh:
        for id_ in id_list:
            fh.write(f"{id_}{os.linesep}")
