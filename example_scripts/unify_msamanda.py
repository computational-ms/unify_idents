from pathlib import Path
from unify_idents.unify import Unify
import sys


def main():
    input_file = (
        Path(__file__).parent.parent / "tests" / "data" / "BSA_msamanda_2_0_0_17442.csv"
    )
    rt_lookup_path = (
        Path(__file__).parent.parent / "tests" / "data" / "BSA_ursgal_lookup.csv.bz2"
    )
    db_path = Path(__file__).parent.parent / "tests" / "data" / "BSA.fasta"

    u = Unify(
        input_file,
        {
            "rt_pickle_name": rt_lookup_path,
            "database": db_path,
            "modifications": [
                "C,fix,any,Carbamidomethyl",
                "M,opt,any,Oxidation",
                "*,opt,Prot-N-term,Acetyl",
            ],
            "Raw file location": "BSA.mzML",
        },
    )

    # u = Unify(
    #     sys.argv[1],
    #     {
    #         "scan_rt_lookup_file": sys.argv[2],
    #         "database": sys.argv[3],
    #         "Modifications": [
    #             "C,fix,any,Carbamidomethyl",
    #             "M,opt,any,Oxidation",
    #             "*,opt,Prot-N-term,Acetyl",
    #         ],
    #         "Raw file location": sys.argv[4],
    #     },
    # )

    df = u.get_dataframe()
    print(df.df)
    df.df.to_csv(
        f'{Path(__file__).parent.parent / "tests" / "data" / "msamanda_unfiy_idents.csv"}'
    )


if __name__ == "__main__":
    main()
