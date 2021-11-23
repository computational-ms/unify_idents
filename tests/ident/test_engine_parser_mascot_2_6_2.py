#!/usr/bin/env python
from pathlib import Path

import pandas as pd
import pytest

from unify_idents.engine_parsers.ident.mascot_2_6_2_parser import (
    Mascot_2_6_2_Parser,
    _get_single_spec_df,
)


def test_engine_parsers_comet_init():
    input_file = pytest._test_path / "data" / "BSA1_mascot_2_6_2.dat"
    rt_lookup_path = pytest._test_path / "data" / "BSA1_ursgal_lookup.csv.bz2"
    db_path = pytest._test_path / "data" / "test_Creinhardtii_target_decoy.fasta"

    parser = Mascot_2_6_2_Parser(
        input_file,
        params={
            "cpus": 2,
            "rt_pickle_name": rt_lookup_path,
            "database": db_path,
            "modifications": [
                {
                    "aa": "M",
                    "type": "opt",
                    "position": "any",
                    "name": "Oxidation",
                },
                {
                    "aa": "C",
                    "type": "fix",
                    "position": "any",
                    "name": "Carbamidomethyl",
                },
                {
                    "aa": "*",
                    "type": "opt",
                    "position": "Prot-N-term",
                    "name": "Acetyl",
                },
            ],
        },
    )


def test_engine_parsers_comet_check_parser_compatibility():
    msgf_parser_class = Mascot_2_6_2_Parser
    input_file = pytest._test_path / "data" / "BSA1_mascot_2_6_2.dat"
    assert msgf_parser_class.check_parser_compatibility(input_file) is True


def test_engine_parsers_comet_check_parser_compatibility_fail_with_omssa_file():
    msgf_parser_class = Mascot_2_6_2_Parser
    input_file = pytest._test_path / "data" / "test_Creinhardtii_QE_pH11_omssa_2_1_9.csv"
    assert msgf_parser_class.check_parser_compatibility(input_file) is False


def test_engine_parsers_comet_check_dataframe_integrity():
    input_file = pytest._test_path / "data" / "BSA1_mascot_2_6_2.dat"
    rt_lookup_path = pytest._test_path / "data" / "BSA1_ursgal_lookup.csv.bz2"
    db_path = pytest._test_path / "data" / "BSA.fasta"

    parser = Mascot_2_6_2_Parser(
        input_file,
        params={
            "cpus": 2,
            "rt_pickle_name": rt_lookup_path,
            "database": db_path,
            "modifications": [
                {
                    "aa": "M",
                    "type": "opt",
                    "position": "any",
                    "name": "Oxidation",
                },
                {
                    "aa": "C",
                    "type": "fix",
                    "position": "any",
                    "name": "Carbamidomethyl",
                },
                {
                    "aa": "*",
                    "type": "opt",
                    "position": "Prot-N-term",
                    "name": "Acetyl",
                },
            ],
        },
    )
    df = parser.unify()

    assert pytest.approx(df["uCalc m/z"].mean()) == 465.30768
    assert pytest.approx(df["Exp m/z"].mean()) == 465.3078

    assert df["Modifications"].str.contains("Acetyl:0").sum() == 3
    assert df["Modifications"].str.contains("Oxidation:").sum() == 0
    assert (
        df["Modifications"].str.count("Carbamidomethyl:")
        == df["Sequence"].str.count("C")
    ).all()
    assert df["Modifications"].str.count(":").sum() == 59
    assert (df["Raw data location"] == "path/for/glory.mzML").all()


def test_get_single_spec_df():
    spec = (
        "query20",
        '"\n\ntitle=BSA1%2e2941%2e2941%2e2\nscans=2941\nrtinseconds=2010.87902832031\nindex=499\ncharge=2+\nmass_min=112.033691\nmass_max=631.522095\nint_min=11.25\nint_max=1.949e+04\nnum_vals=133\nnum_used1=-1\nIons1=129.049026:3748,244.139252:1.387e+04,315.724304:1.949e+04,470.263153:1.047e+04,515.171631:1.081e+04,630.382629:1.898e+04,133.205566:887.7,289.091186:4416,402.107849:2390,487.249298:3452,516.346985:1112,631.522095:2070,147.092804:828.9,226.108429:3390,351.299927:946.7,424.381775:807.6,595.231201:107.8,612.390442:1626,130.232208:561.5,306.911011:2655,387.218658:621.3,497.235016:724.5,594.108459:100.7,613.442871:450,209.203354:496.4,311.232941:949.7,339.289429:602.8,471.278625:656.2,584.450867:59.95,614.067810:24.01,183.196442:348.1,246.180817:825.5,353.655365:267.9,413.121246:253.7,543.242188:45.01,116.150803:307.9,245.091186:774.2,331.369751:241.6,498.236542:248.7,539.399719:39.76,146.116867:275.9,274.044250:558,406.513489:238.2,488.305695:171.2,211.179108:185.6,227.054504:442.8,325.410675:197.2,499.225311:133.3,198.115265:124.9,260.864166:348.2,352.833832:175.4,502.027008:126.6,112.033691:23.2,116.836945:14.23,126.112038:38.97,135.970749:26.17,150.372131:35.92,155.268936:106.4,159.086746:37.79,160.938324:35.79,169.198685:29.97,178.182220:45.89,180.375137:71.85,181.199539:39.61,196.922821:101.2,199.325928:63.47,202.350220:58.42,204.328705:95.82,210.025406:38.9,214.150085:62.69,214.992752:80.54,228.156738:20.31,229.194916:134.2,230.055725:73.07,232.120666:64.12,236.289551:59.01,239.429352:203.7,240.094513:20.23,248.428223:21.72,250.192108:35.81,257.314026:100.4,259.325653:50.72,263.221802:14.14,266.196289:66.5,268.133240:61.32,271.065674:100.8,272.133637:32.16,275.988037:323,283.323547:65.88,284.492401:18.82,286.912445:66.83,290.186676:44.81,291.217407:68.01,293.195190:24.01,296.132202:41.13,298.968048:105,307.886475:209.2,309.374084:20.14,312.279083:105.7,317.097473:120.2,317.919006:34.97,318.968933:123,319.995544:48.06,321.120331:63.42,322.493591:85.83,323.225281:38.16,327.023773:92.02,329.155090:84.8,340.023987:40.51,342.130951:26.17,343.166290:31.97,345.129211:82.29,355.029694:42,384.251739:39.83,385.280487:31.47,386.294708:101.6,388.347900:116,389.094025:72.61,403.279999:81.47,407.245544:68.93,425.129242:24.47,426.360870:41.17,427.261322:105.3,428.050812:17.16,430.350250:80.5,437.316315:82.45,444.145294:55.92,456.027191:40.51,458.869354:18.73,477.371155:11.25,480.504303:94.21,484.294159:119.6,486.008575:61.18\n--gc0p4Jq0M2Yt08jU534c0p\n',
        [
            '0,757.415634,-0.000498,5,GACLLPK,18,000000000,35.26,0002001000000000000,0,0;"sp|P02769|ALBU_BOVIN":0:198:204:1;subst;None'
        ],
    )
    ref_dict = {
        "Exp m/z": None,
        "Calc m/z": None,
        "Spectrum Title": None,
        "Raw data location": "/Users/av568207/dev/ursgal2/example_data/BSA1.mgf",
        "Search Engine": "mascot_2_6_2",
        "Spectrum ID": None,
        "Modifications": None,
        "Retention Time (s)": None,
        "Mascot:Score": None,
    }

    result = _get_single_spec_df(ref_dict, spec)
    assert isinstance(result, pd.DataFrame)
    assert (
        result.values
        == [
            [
                "757.415634",
                None,
                "BSA1%2e2941%2e2941%2e2",
                "/Users/av568207/dev/ursgal2/example_data/BSA1.mgf",
                "mascot_2_6_2",
                "2941",
                "000000000",
                None,
                "35.26",
                "2",
                "5",
                "GACLLPK",
                "None",
            ]
        ]
    ).all()


def test_translate_opt_mods():
    input_file = pytest._test_path / "data" / "BSA1_mascot_2_6_2.dat"
    rt_lookup_path = pytest._test_path / "data" / "BSA1_ursgal_lookup.csv.bz2"
    db_path = pytest._test_path / "data" / "BSA.fasta"

    parser = Mascot_2_6_2_Parser(
        input_file,
        params={
            "cpus": 2,
            "rt_pickle_name": rt_lookup_path,
            "database": db_path,
            "modifications": [
                {
                    "aa": "M",
                    "type": "opt",
                    "position": "any",
                    "name": "Oxidation",
                },
                {
                    "aa": "C",
                    "type": "fix",
                    "position": "any",
                    "name": "Carbamidomethyl",
                },
                {
                    "aa": "*",
                    "type": "opt",
                    "position": "Prot-N-term",
                    "name": "Acetyl",
                },
            ],
        },
    )
    raw_mod = "2000000000"
    converted = parser._translate_opt_mods(raw_mod)
    assert converted == ";Acetyl:0;"
