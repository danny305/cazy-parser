from pathlib import Path
import argparse
import logging
import sys
import time

from cazy_parser import ENZYME_LIST
from cazy_parser.modules.fasta import dump_fastas, dump_id_list
from cazy_parser.modules.html import retrieve_cazy_metadata
from cazy_parser.modules.metadata import dump_metadata
from cazy_parser.version import VERSION

log = logging.getLogger("cazylog")
ch = logging.StreamHandler()
formatter = logging.Formatter(" [%(asctime)s %(lineno)d %(levelname)s] %(message)s")
ch.setFormatter(formatter)
log.addHandler(ch)
log.setLevel("DEBUG")


# ====================================================================================#
# Main code
def main():
    """Main function."""

    ap = argparse.ArgumentParser()

    group = ap.add_mutually_exclusive_group(required=True)
    group.add_argument("--fasta", action="store_true")
    group.add_argument("--metadata", action="store_true")

    ap.add_argument(
        "enzyme_class",
        choices=["GH", "GT", "PL", "CA", "AA"],
    )

    ap.add_argument("-f", "--family", type=int, default=None)

    ap.add_argument("-s", "--subfamily", type=int, default=None)

    ap.add_argument("-c", "--characterized", action="store_true", default=False)

    ap.add_argument("-d", "--out-dir", type=Path, default=Path.cwd() / "fasta")

    ap.add_argument("-i", "--genbank-id", action="store_true")

    ap.add_argument("--debug-level", choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"], default="INFO")


    ap.add_argument(
        "-v",
        "--version",
        help="show version",
        action="version",
        version=f"Running {ap.prog} v{VERSION}",
    )

    args = ap.parse_args()

    assert args.debug_level in ["DEBUG", "INFO", "WARNING", "ERROR" "CRITICAL"], f"Invalid debug level {args.debug_level}"

    log.setLevel(args.debug_level)

    if args.enzyme_class not in ENZYME_LIST:
        logging.error(f"Enzyme class {args.enzyme_class} not supported")
        sys.exit()
    else:
        enzyme_name = ENZYME_LIST[args.enzyme_class]

    log.info("-" * 42)
    log.info("")
    log.info("┌─┐┌─┐┌─┐┬ ┬   ┌─┐┌─┐┬─┐┌─┐┌─┐┬─┐")
    log.info("│  ├─┤┌─┘└┬┘───├─┘├─┤├┬┘└─┐├┤ ├┬┘")
    log.info(f"└─┘┴ ┴└─┘ ┴    ┴  ┴ ┴┴└─└─┘└─┘┴└─ v{VERSION}")
    log.info("")
    log.info("-" * 42)

    # id_list = retrieve_genbank_ids(
    #     enzyme_name, args.family, args.subfamily, args.characterized
    # )

    cazy_metdata = retrieve_cazy_metadata(
        enzyme_name, 
        args.family, 
        args.subfamily, 
        args.characterized,
        "genbank"
    )

    id_list = [element["genbank"] for element in cazy_metdata]

    output_fname = f"{args.enzyme_class}"
    if args.family:
        output_fname += f"{args.family}"
    if args.subfamily:
        output_fname += f"_{args.subfamily}"
    if args.characterized:
        output_fname += "_characterized"

    today = time.strftime("%m%d%y")
    output_fname += f"_{today}"

    gb_dir = args.out_dir / "genbank" 
    gb_dir.mkdir(0o774, parents=True, exist_ok=True)

    fasta_dir = args.out_dir / "fasta"
    meta_dir = args.out_dir / "metadata"

    gb_id_fname = gb_dir / f"{output_fname}.txt"
    if args.fasta:
        try:
            fasta_dir.mkdir(0o774, parents=True, exist_ok=True) 
            fasta_out = fasta_dir / f"{output_fname}.fasta"
            dump_fastas(id_list, fasta_out)

        except Exception as e:
            log.debug(e)
            log.warning(
                "Could not fetch the fasta sequences, dumping the sequence IDs instead."
            )
            log.warning(
                "This is probably due to the NCBI server being inaccessible. Please try again later or manually download the sequences from NCBI"
            )
            log.warning(
                f"Please upload {gb_id_fname} to `https://www.ncbi.nlm.nih.gov/sites/batchentrez` to download the sequences"
            )
        finally: 
            dump_id_list(id_list, gb_id_fname)

    else:
        meta_dir.mkdir(0o774, parents=True, exist_ok=True)
        meta_out = meta_dir / f"{output_fname}.csv"
        dump_metadata(cazy_metdata, meta_out)


if __name__ == "__main__":
    sys.exit(main())
