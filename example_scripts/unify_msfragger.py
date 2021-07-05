from pathlib import Path
from unify_idents.unify import Unify
import sys


def main(input=None):

    u = Unify(
        sys.argv[1],
        {
            "scan_rt_lookup_file": sys.argv[2],
            "database": sys.argv[3],
            "Modifications": [
                "C,fix,any,Carbamidomethyl",
                "M,opt,any,Oxidation",
                "*,opt,Prot-N-term,Acetyl",
            ],
            "Raw file location": sys.argv[4],
        },
    )

    df = u.get_dataframe()
    print(df.df)
    df.df.to_csv("/Users/cellzome/Data/StSchulze/ursgal_1_vs_2/unfiy_idents.csv")


if __name__ == "__main__":
    main()
