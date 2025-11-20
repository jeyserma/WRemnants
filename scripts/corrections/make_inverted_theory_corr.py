import pickle

import lz4.frame

from utilities import common, parsing
from wums import logging, output_tools


def parse_args():
    parser = parsing.base_parser()
    parser.add_argument(
        "-c",
        "--corrName",
        type=str,
        help="Reference corr (usually a correction to another prediction+unc.)",
        required=True,
    )
    parser.add_argument(
        "-f",
        "--corrFile",
        type=str,
        help="Reference corr (usually a correction to another prediction+unc.)",
        required=True,
    )
    parser.add_argument(
        "-n",
        "--newCorrName",
        type=str,
        help="Name of the new correction",
        required=True,
    )
    parser.add_argument(
        "--outpath",
        type=str,
        default=f"{common.data_dir}/TheoryCorrections",
        help="Output path",
    )

    return parser.parse_args()


def main():
    args = parse_args()
    logger = logging.setup_logger(__file__, args.verbose, args.noColorLogger)

    ref = pickle.load(lz4.frame.open(args.corrFile))

    proc = "Z" if "CorrZ" in args.corrFile else "W"

    refcorr = ref[proc][args.corrName]

    vals = refcorr.values()
    refcorr[...] = 1.0 / (vals + (vals == 0))

    outname = args.newCorrName if "dataPtll" not in args.newCorrName else args.corrName
    output_dict = {
        args.newCorrName: refcorr,
        **ref[proc],
    }
    meta_dict = {
        "ref_file": ref["file_meta_data"],
    }
    outfile = f"{args.outpath}/{args.newCorrName}Corr{proc}.pkl.lz4"

    output_tools.write_lz4_pkl_output(
        outfile, proc, output_dict, common.base_dir, args, meta_dict
    )


if __name__ == "__main__":
    main()
