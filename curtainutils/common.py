from typing import List

import requests
import pandas as pd
import json
from copy import deepcopy
from uniprotparser.betaparser import UniprotSequence

curtain_base_de_form = {
    "_reverseFoldChange": False,
    "_comparisonSelect": [],
    "_comparison": "",
    "_primaryIDs": "",
    "_geneNames": "",
    "_foldChange": "",
    "_transformFC": False,
    "_significant": "",
    "_transformSignificant": False,
}

curtain_base_raw_form = {
    "_primaryIDs": "",
    "_samples": [],
}

curtain_base_input_file = {
    "df": "",
    "filename": "",
    "other": {},
    "originalFile": ""
}

curtain_base_project_form = {
    "title": "",
    "projectDescription": "",
    "organisms": [{"name": ""}],
    "organismParts": [{"name": ""}],
    "cellTypes": [{"name": ""}],
    "diseases": [{"name": ""}],
    "sampleProcessingProtocol": "",
    "dataProcessingProtocol": "",
    "identifiedPTMStrings": [{"name": ""}],
    "instruments": [{"name": "", "cvLabel": "MS"}],
    "msMethods": [{"name": ""}],
    "projectTags": [{"name": ""}],
    "quantificationMethods": [{"name": ""}],
    "species": [{"name": ""}],
    "sampleAnnotations": {},
    "_links": {"datasetFtpUrl": {"href": ""}, "files": {"href": ""}, "self": {"href": ""}},
    "affiliations": [{"name": ""}],
    "hasLink": False,
    "authors": [],
    "accession": "",
    "softwares": [{"name": ""}],
}



curtain_base_settings = {
    "fetchUniprot": True,
    "sampleMap": {},
    "barchartColorMap": {},
    "pCutoff": 0.05,
    "log2FCCutoff": 0.6,
    "description": "",
    "uniprot": False,
    "colorMap": {},
    "backGroundColorGrey": False,
    "selectedComparison": [],
    "version": 2,
    "currentID": "",
    "fdrCurveText": "",
    "fdrCurveTextEnable": False,
    "prideAccession": "",
    "project": curtain_base_project_form,
    "sampleOrder": {},
    "sampleVisible": {},
    "conditionOrder": [],
    "volcanoAxis": {
        "minX": None,
        "maxX": None,
        "minY": None,
        "maxY": None
    },
    "textAnnotation": {},
    "volcanoPlotTitle": "",
    "visible": {},
    "defaultColorList": [
        "#fd7f6f",
        "#7eb0d5",
        "#b2e061",
        "#bd7ebe",
        "#ffb55a",
        "#ffee65",
        "#beb9db",
        "#fdcce5",
        "#8bd3c7",
    ],
    "scatterPlotMarkerSize": 10,
    "rankPlotColorMap": {},
    "rankPlotAnnotation": {},
    "legendStatus": {},
    "stringDBColorMap": {
        "Increase": "#8d0606",
        "Decrease": "#4f78a4",
        "In dataset": "#ce8080",
        "Not in dataset": "#676666"
    },
    "interactomeAtlasColorMap": {
        "Increase": "#a12323",
        "Decrease": "#16458c",
        "HI-Union": "rgba(82,110,194,0.96)",
        "Literature": "rgba(181,151,222,0.96)",
        "HI-Union and Literature": "rgba(222,178,151,0.96)",
        "Not found": "rgba(25,128,128,0.96)",
        "No change": "rgba(47,39,40,0.96)",
    },
    "proteomicsDBColor": "#ff7f0e",
    "networkInteractionSettings": {
        "Increase": "rgba(220,169,0,0.96)",
        "Decrease": "rgba(220,0,59,0.96)",
        "StringDB": "rgb(206,128,128)",
        "No change": "rgba(47,39,40,0.96)",
        "Not significant": "rgba(255,255,255,0.96)",
        "Significant": "rgba(252,107,220,0.96)",
        "InteractomeAtlas": "rgb(73,73,101)",
    },
    "plotFontFamily": "Arial",
    "networkInteractionData": [],
    "enrichrGeneRankMap": {},
    "enrichrRunList": [],
    "customVolcanoTextCol": ""
}

curtain_base_payload = {
    "raw": "",
    "rawForm": curtain_base_raw_form,
    "differentialForm": curtain_base_de_form,
    "processed": "",
    "password": "",
    "selections": [],
    "selectionsMap": {},
    "selectionsName": [],
    "settings": curtain_base_settings,
    "fetchUniprot": True,
    "annotatedData": {},
}
def read_fasta(fasta_file: str) -> pd.DataFrame:
    fasta_dict = {}
    with open(fasta_file, 'r') as f:
        current_acc = ""
        count = 0
        for line in f:
            if line.startswith('>'):
                count += 1
                acc = UniprotSequence(line.strip(), True)

                if acc.accession:
                    fasta_dict[str(acc)] = ""
                    current_acc = str(acc)
                else:
                    fasta_dict[line.strip().replace(">", "")] = ""
                    current_acc = line.strip().replace(">", "")

            else:
                fasta_dict[current_acc] += line.strip()


    return pd.DataFrame([[k, fasta_dict[k]] for k in fasta_dict], columns=["Entry", "Sequence"])


def create_curtain_session_payload(
        de_file: str,
        raw_file: str,
        fc_col: str,
        transform_fc: bool,
        transform_significant: bool,
        reverse_fc: bool,
        p_col: str,
        comp_col: str,
        comp_select: List[str],
        primary_id_de_col: str,
        primary_id_raw_col:str,
        sample_cols: List[str], **kwargs) -> dict:
    payload = deepcopy(curtain_base_payload)
    with open(de_file, "rt") as f, open(raw_file, "rt") as f2:
        payload["processed"] = f.read()
        payload["raw"] = f2.read()

    payload["differentialForm"]["_foldChange"] = fc_col
    payload["differentialForm"]["_significant"] = p_col
    payload["differentialForm"]["_comparison"] = comp_col
    payload["differentialForm"]["_comparisonSelect"] = comp_select
    payload["differentialForm"]["_primaryIDs"] = primary_id_de_col
    payload["rawForm"]["_primaryIDs"] = primary_id_raw_col
    payload["rawForm"]["_samples"] = sample_cols

    assert type(transform_fc) == bool
    assert type(transform_significant) == bool
    assert type(reverse_fc) == bool

    payload["differentialForm"]["_transformFC"] = transform_fc
    payload["differentialForm"]["_transformSignificant"] = transform_significant
    payload["differentialForm"]["_reverseFoldChange"] = reverse_fc

    if payload["differentialForm"]["_comparison"] == "":
        payload["differentialForm"]["_comparison"] = "CurtainSetComparison"

    if len(payload["differentialForm"]["_comparisonSelect"]) == 0:
        payload["differentialForm"]["_comparisonSelect"] = ["1"]

    assert len(sample_cols) > 0
    conditions = []
    color_position = 0
    sample_map = {}
    color_map = {}
    for i in sample_cols:
        name_array = i.split(".")
        replicate = name_array[-1]
        condition = ".".join(name_array[:-1])
        if condition not in conditions:
            conditions.append(condition)
            if color_position >= len(payload["settings"]["defaultColorList"]):
                color_position = 0
            color_map[condition] = payload["settings"]["defaultColorList"][color_position]
            color_position += 1
        if condition not in payload["settings"]["sampleOrder"]:
            payload["settings"]["sampleOrder"][condition] = []
        if i not in payload["settings"]["sampleOrder"][condition]:
            payload["settings"]["sampleOrder"][condition].append(i)
        if i not in payload["settings"]["sampleVisible"]:
            payload["settings"]["sampleVisible"][i] = True

        sample_map[i] = {"condition": condition, "replicate": replicate, "name": i}
    payload["settings"]["sampleMap"] = sample_map
    payload["settings"]["colorMap"] = color_map
    payload["settings"]["conditionOrder"] = conditions

    return payload


def post_curtain_session(payload: dict, file: dict, url: str, **kwargs) -> dict:
    #headers = {'Content-Type': 'multipart/form-data', 'Accept': 'application/json'}
    file = {'file': ('curtain-settings.json', json.dumps(file))}
    r = requests.post(url, data=payload, files=file)
    return r
