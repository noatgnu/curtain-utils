import os.path

import pandas as pd
from uniprotparser.betaparser import UniprotSequence
import base64
from cryptography.hazmat.backends import default_backend
from cryptography.hazmat.primitives.asymmetric import rsa
from cryptography.hazmat.primitives.ciphers import Cipher, algorithms, modes
from cryptography.hazmat.primitives import serialization
from cryptography.hazmat.primitives import hashes
from cryptography.hazmat.primitives.asymmetric import padding
from cryptography.hazmat.primitives.padding import PKCS7

curtain_ptm_de_form = {
    "_reverseFoldChange": False,
    "_comparisonSelect": "",
    "_comparison": "",
    "_primaryIDs": "",
    "_geneNames": "",
    "_foldChange": "",
    "_transformFC": False,
    "_significant": "",
    "_transformSignificant": False,
    "_accession": "",
    "_position": "",
    "_positionPeptide": "",
    "_peptideSequence": "",
    "_score": "",
}

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

curtain_base_input_file = {"df": "", "filename": "", "other": {}, "originalFile": ""}

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
    "_links": {
        "datasetFtpUrl": {"href": ""},
        "files": {"href": ""},
        "self": {"href": ""},
    },
    "affiliations": [{"name": ""}],
    "hasLink": False,
    "authors": [],
    "accession": "",
    "softwares": [{"name": ""}],
}


curtain_base_settings = {
    "fetchUniProt": True,
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
    "volcanoAxis": {"minX": None, "maxX": None, "minY": None, "maxY": None},
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
        "Not in dataset": "#676666",
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
    "customVolcanoTextCol": "",
}

curtain_ptm_settings = {
    "fetchUniProt": True,
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
    "volcanoAxis": {"minX": None, "maxX": None, "minY": None, "maxY": None},
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
        "Not in dataset": "#676666",
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
    "customVolcanoTextCol": "",
    "variantCorrection": {},
    "customSequences": {},
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
    "fetchUniProt": True,
    "annotatedData": {},
}


def read_fasta(fasta_file: str) -> pd.DataFrame:
    fasta_dict = {}
    with open(fasta_file, "r") as f:
        current_acc = ""
        count = 0
        for line in f:
            if line.startswith(">"):
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

    return pd.DataFrame(
        [[k, fasta_dict[k]] for k in fasta_dict], columns=["Entry", "Sequence"]
    )


def create_ptm_db(
    input_file: str,
    uniprot_acc_col: str,
    peptide_seq_col: str = "",
    peptide_pos_col: str = "",
    modified_residue_is_lower_case: bool = False,
    output_file: str = "",
    peptide_start_col: str = "",
    parse_fasta: bool = False,
    fasta_file: str = "",
):
    df = pd.read_csv(input_file, sep="\t")

    if modified_residue_is_lower_case:
        # parse position of modified residue from peptide sequence column by detecting lower case letter
        df["Position"] = df[peptide_seq_col].apply(
            lambda x: x.index([i for i in x if i.islower()][0]) + 1
        )
        if peptide_start_col != "":
            df["Position"] += df[peptide_start_col] - 1
        else:
            parse_fasta = True

    return df


def generate_RSA_key_pair(folder_path: str = "."):
    private_key = rsa.generate_private_key(
        public_exponent=65537, key_size=2048, backend=default_backend()
    )
    public_key = private_key.public_key()
    private_pem = private_key.private_bytes(
        encoding=serialization.Encoding.PEM,
        format=serialization.PrivateFormat.PKCS8,
        encryption_algorithm=serialization.NoEncryption(),
    )

    with open(os.path.join(folder_path, "private_key.pem"), "wb") as f:
        f.write(private_pem)

    public_pem = public_key.public_bytes(
        encoding=serialization.Encoding.PEM,
        format=serialization.PublicFormat.SubjectPublicKeyInfo,
    )
    with open(os.path.join(folder_path, "public_key.pem"), "wb") as f:
        f.write(public_pem)


# a function to generate a aes key and encrypt it with an RSA public key
def encrypt_AES_key(aes_key: bytes, public_key: str) -> bytes:
    public_key = serialization.load_pem_public_key(
        public_key, backend=default_backend()
    )
    encrypted = public_key.encrypt(
        aes_key,
        padding.OAEP(
            mgf=padding.MGF1(algorithm=hashes.SHA256()),
            algorithm=hashes.SHA256(),
            label=None,
        ),
    )
    return encrypted


# a function to encrypt a large string with an AES key in CTR mode and return the encrypted string and the initialization vector
def encrypt_AES_string(aes_key: bytes, string: str) -> (bytes, bytes):
    iv = os.urandom(16)
    cipher = Cipher(algorithms.AES(aes_key), modes.CTR(iv), backend=default_backend())
    encryptor = cipher.encryptor()
    ct = encryptor.update(string.encode("utf-8")) + encryptor.finalize()
    return ct, iv
