from pathlib import Path
from typing import Literal
from time import time
import logging
import os

from Bio import Entrez, SeqIO
from tqdm import tqdm

Entrez.email = "your-email-here@example.org"

log = logging.getLogger("cazylog")

from cazy_parser.modules.uniprot import fetch_uniprot_sequence


def flatten(data: list[str]|list[list[str]]) -> list[str]:
    flat_list = []
    for item in data:
        if isinstance(item, list):
            flat_list.extend(flatten(item))
        elif ":" in item:
            flat_list.extend(item.split(":"))
        else:
            flat_list.append(item)
    return flat_list


def download_genbank_fastas(id_list: list[str]) -> dict[str:str]:
    """
    Download fasta files from NCBI.

    Parameters
    ----------
    id_list : dictionary
        List of genbank ids.

    Returns
    -------
    fasta_list : dictionary
        dictionary with id as key and fasta string as value.

    """
    id_list = {item for item in id_list if item.replace("\xa0","")}
    start_t = time()
    log.info(f"Dowloading {len(id_list)} fasta sequences...")
    fastas = {}
    handle = Entrez.efetch(db="protein", id=id_list, rettype="fasta", retmode="text")
    log.info(f"Entrez api call: {handle.url}")
    for seq_record in tqdm(SeqIO.parse(handle, "fasta"), total=len(id_list)):
        fasta_str = (
            f">{seq_record.description}{os.linesep}{seq_record.seq}{os.linesep*2}"
        )
        fastas[seq_record.id] = fasta_str

    log.info(f"Downloaded {len(fastas.keys())} fasta sequences in {time()-start_t:.2f} seconds")

    if len(fastas) != len(id_list):
        log.warning(
            f"Only {len(fastas)} fasta sequences were downloaded out of {len(id_list)}" 
        )
        
    return fastas


def dump_fastas(id_list: list[str], id_type: Literal["genbank", "uniprot"], output_f: Path):
    """
    Save the fasta strings to a file.

    Parameters
    ----------
    id_list : list
        List of genbank ids.
    output_f : str
        Path to the output file.

    """
    assert id_type in {"genbank", "uniprot"}, f"{id_type} is not a valid id type"
    assert isinstance(output_f, Path), f"{output_f} is of type {type(output_f)}"

    fastas = download_genbank_fastas(id_list) if id_type == "genbank" else download_uniprot_fastas(id_list)
    log.info(f"Dumping fasta sequences to file {output_f.resolve()}")
    with output_f.open("w") as fh:
        for fasta in fastas.values():
            fh.write(fasta)

    log.info(f"IDs missing a fasta: {set(id_list) - set(fastas.keys())}")


def dump_id_list(id_list: list[str], output_f: Path) -> None:
    """Save the id list to a file."""
    assert isinstance(output_f, Path), f"{output_f} is of type {type(output_f)}"
    log.info(f"Dumping ids to file {output_f.resolve()}")
    with output_f.open("w") as fh:
        for id_ in id_list:
            fh.write(f"{id_}{os.linesep}")


def download_uniprot_fastas(uniprot_ids: list[str]) -> dict[str:str]:
    """
    Download fasta files from UNIPROT.

    Parameters
    ----------
    id_list : list
        List of genbank ids.

    Returns
    -------
    fasta_list : dictionary
        dictionary with id as key and fasta string as value.
    """

    uniprot_ids = {item for item in uniprot_ids if item.replace("\xa0","")}

    fastas = {}
    for uniprot_id in tqdm(uniprot_ids):
        header, sequence = fetch_uniprot_sequence(uniprot_id)
        if header and sequence:
            fasta_str = f"{header}{os.linesep}{sequence}{os.linesep*2}"
            fastas[uniprot_id] = fasta_str
        else:
            log.warning(f"Failed to fetch {uniprot_id}")

    return fastas