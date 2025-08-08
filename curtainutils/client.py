import io
import json
import re
from copy import deepcopy
from typing import List, Dict, Optional, Generator
import time
import pandas as pd
import requests
import seaborn as sns
from matplotlib import pyplot as plt

from curtainutils.common import curtain_base_payload
from uniprotparser.betaparser import UniprotParser, UniprotSequence
import numpy as np


class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super(NumpyEncoder, self).default(obj)


class CurtainUniprotData:
    def __init__(self, uniprot: Dict):
        """
        Initialize CurtainUniprotData with UniProt data.

        Args:
            uniprot: Dictionary containing UniProt data
        """
        self.accMap = {i[0]: i[1] for i in uniprot["accMap"]["value"]}
        self.dataMap = {i[0]: i[1] for i in uniprot["dataMap"]["value"]}
        self.db = self._db_to_df(uniprot["db"])
        self.organism = uniprot.get("organism")
        self.geneNameToAcc = uniprot.get("geneNameToAcc", {})
        self.geneNameToPrimary = uniprot.get("geneNameToPrimary", {})
        self.results = uniprot["results"]

    def _db_to_df(self, db: Dict) -> pd.DataFrame:
        """Convert database dictionary to DataFrame"""
        data = [i[1] for i in db["value"]]
        return pd.DataFrame(data)

    def get_uniprot_data_from_pi(self, primary_id: str) -> Optional[pd.Series]:
        """
        Get UniProt data for a primary ID.

        Args:
            primary_id: Primary identifier to search for

        Returns:
            DataFrame row containing the UniProt data or None if not found
        """
        return CurtainUniprotData.get_uniprot_data_from_pi_sta(
            primary_id, self.accMap, self.dataMap, self.db
        )

    @staticmethod
    def get_uniprot_data_from_pi_sta(
        primary_id: str, accMap: Dict, dataMap: Dict, db: pd.DataFrame
    ) -> Optional[pd.Series]:
        if primary_id in accMap:
            acc_match_list = accMap[primary_id]
            if acc_match_list:
                if isinstance(acc_match_list, str):
                    acc_match_list = [acc_match_list]
                for acc in acc_match_list:
                    if acc in dataMap:
                        ac = dataMap[acc]
                        filter_db = db[db["From"] == ac]
                        if not filter_db.empty:
                            return filter_db.iloc[0]
        return None


class CurtainClient:
    def __init__(self, base_url: str, api_key: str = ""):
        """
        Initialize Curtain client for API interaction.

        Args:
            base_url: Base URL for Curtain API
            api_key: Optional API key for authentication
        """
        self.base_url = base_url.rstrip("/")
        self.request_session = requests.Session()
        self.refresh_token = ""
        self.access_token = ""
        self.api_key = api_key

    def _get_auth_headers(self) -> Dict[str, str]:
        """Get authentication headers if API key is provided"""
        return {"X-Api-Key": self.api_key} if self.api_key else {}

    def get_data_filter_list(self) -> Generator[List, None, None]:
        """
        Get data filter list from API with pagination support.

        Yields:
            Results from each page of the API response

        Raises:
            requests.HTTPError: If API request fails
        """
        headers = self._get_auth_headers()
        r = requests.get(f"{self.base_url}/data_filter_list/", headers=headers)
        r.raise_for_status()

        res = r.json()
        yield res["results"]

        while res.get("next"):
            r = requests.get(res["next"], headers=headers)
            r.raise_for_status()
            res = r.json()
            yield res["results"]

    def post_curtain_session(self, payload: Dict, file: Dict) -> str:
        """
        Post a new Curtain session.

        Args:
            payload: Session metadata
            file: Session file data

        Returns:
            Link ID of the created session

        Raises:
            requests.HTTPError: If API request fails
        """
        file_data = {
            "file": ("curtain-settings.json", json.dumps(file, cls=NumpyEncoder))
        }
        headers = self._get_auth_headers()

        r = requests.post(
            f"{self.base_url}/curtain/", data=payload, files=file_data, headers=headers
        )
        r.raise_for_status()

        return r.json()["link_id"]

    def curtain_stats_summary(self, last_n_days: int = 30) -> None:
        """
        Display stats summary for Curtain sessions.

        Args:
            last_n_days: Number of days to include in the summary
        """
        req = requests.get(f"{self.base_url}/stats/summary/{last_n_days}/")
        req.raise_for_status()
        data = req.json()

        session_download_per_day = pd.DataFrame(data["session_download_per_day"])
        session_download_per_day.sort_values(by="date", inplace=True)

        session_created_per_day = pd.DataFrame(data["session_created_per_day"])
        session_created_per_day.sort_values(by="date", inplace=True)

        fig, ax = plt.subplots(2, 1, figsize=(20, 10))
        sns.barplot(data=session_download_per_day, x="date", y="downloads", ax=ax[0])
        sns.barplot(data=session_created_per_day, x="date", y="count", ax=ax[1])

        ax[0].set_title("Curtain Session Download per day")
        ax[1].set_title("Curtain Session Creation per day")

        for axis in fig.axes:
            plt.sca(axis)
            plt.xticks(rotation=90)

        fig.tight_layout()
        plt.show()

    def download_curtain_session(
        self,
        link_id: str,
        retries: int = 3,
        show_progress: bool = False,
        progress_callback=None,
    ) -> Optional[Dict]:
        """
        Download Curtain session data with retry, optional progress display, and callback support.

        Args:
            link_id: ID of the session to download
            retries: Number of times to retry on failure (default: 3)
            show_progress: If True, display download progress bar (default: False)
            progress_callback: Optional callback function that receives (downloaded_bytes, total_bytes, percentage)

        Returns:
            Session data or None if download fails
        """

        link = f"{self.base_url}/curtain/{link_id}/download/token=/"

        for attempt in range(retries + 1):
            try:
                req = requests.get(link)

                if req.status_code == 200:
                    data = req.json()
                    if "url" in data:
                        # Download the actual data with optional progress
                        for dl_attempt in range(retries + 1):
                            try:
                                if show_progress:
                                    try:
                                        from tqdm import tqdm
                                    except ImportError:
                                        print(
                                            "tqdm not available, downloading without progress bar..."
                                        )
                                        show_progress = False

                                if show_progress or progress_callback:
                                    with requests.get(
                                        data["url"], stream=True
                                    ) as result:
                                        if result.status_code == 200:
                                            total = int(
                                                result.headers.get("content-length", 0)
                                            )
                                            downloaded = 0
                                            chunks = []

                                            # Initialize progress bar if needed
                                            pbar = None
                                            if show_progress:
                                                pbar = tqdm(
                                                    total=total,
                                                    unit="B",
                                                    unit_scale=True,
                                                    desc="Downloading",
                                                )

                                            try:
                                                for chunk in result.iter_content(
                                                    chunk_size=8192
                                                ):
                                                    if chunk:
                                                        chunks.append(chunk)
                                                        downloaded += len(chunk)

                                                        # Update progress bar
                                                        if pbar:
                                                            pbar.update(len(chunk))

                                                        # Call progress callback if provided
                                                        if progress_callback:
                                                            percentage = (
                                                                (
                                                                    downloaded
                                                                    / total
                                                                    * 100
                                                                )
                                                                if total > 0
                                                                else 0
                                                            )
                                                            progress_callback(
                                                                downloaded,
                                                                total,
                                                                percentage,
                                                            )
                                            finally:
                                                if pbar:
                                                    pbar.close()

                                            content = b"".join(chunks)
                                            return json.loads(content.decode())
                                        elif dl_attempt < retries:
                                            time.sleep(
                                                2**dl_attempt
                                            )  # Exponential backoff
                                            continue
                                        else:
                                            return None
                                else:
                                    result = requests.get(data["url"])
                                    if result.status_code == 200:
                                        return result.json()
                                    elif dl_attempt < retries:
                                        time.sleep(2**dl_attempt)  # Exponential backoff
                                        continue
                                    else:
                                        return None

                            except (
                                requests.exceptions.RequestException,
                                json.JSONDecodeError,
                            ) as e:
                                if dl_attempt < retries:
                                    time.sleep(2**dl_attempt)  # Exponential backoff
                                    continue
                                return None
                    else:
                        return data
                elif attempt < retries:
                    time.sleep(2**attempt)  # Exponential backoff
                    continue
                else:
                    return None

            except requests.exceptions.RequestException as e:
                if attempt < retries:
                    time.sleep(2**attempt)  # Exponential backoff
                    continue
                return None

        return None

    def download_sessions_list(self, url_list: List[str]) -> None:
        """
        Download multiple Curtain sessions and save to files.

        Args:
            url_list: List of session IDs to download
        """
        for session_id in url_list:
            result = self.download_curtain_session(session_id)
            if result:
                with open(f"{session_id}_de.txt", "wt") as f:
                    f.write(result["processed"])
                with open(f"{session_id}_raw.txt", "wt") as f:
                    f.write(result["raw"])

    def _prepare_common_payload(
        self,
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
        primary_id_raw_col: str,
        sample_cols: List[str],
        description: str = "",
    ) -> Dict:
        """
        Prepare common payload elements for both standard and PTM Curtain sessions.

        Returns:
            Base payload with common elements configured
        """
        payload = deepcopy(curtain_base_payload)

        # Read files
        with open(de_file, "rt") as f, open(raw_file, "rt") as f2:
            payload["processed"] = f.read()
            payload["raw"] = f2.read()
        payload["settings"]["description"] = description
        # Set differential form settings
        payload["differentialForm"]["_foldChange"] = fc_col
        payload["differentialForm"]["_significant"] = p_col
        payload["differentialForm"]["_comparison"] = comp_col or "CurtainSetComparison"
        payload["differentialForm"]["_comparisonSelect"] = comp_select or ["1"]
        payload["differentialForm"]["_primaryIDs"] = primary_id_de_col
        payload["rawForm"]["_primaryIDs"] = primary_id_raw_col
        payload["rawForm"]["_samples"] = sample_cols

        # Set transformation flags
        payload["differentialForm"]["_transformFC"] = transform_fc
        payload["differentialForm"]["_transformSignificant"] = transform_significant
        payload["differentialForm"]["_reverseFoldChange"] = reverse_fc

        # Configure sample and color settings
        self._configure_sample_settings(payload, sample_cols)

        return payload

    def _configure_sample_settings(self, payload: Dict, sample_cols: List[str]) -> None:
        """
        Configure sample-related settings in the payload.

        Args:
            payload: Payload to modify
            sample_cols: List of sample column names
        """
        assert len(sample_cols) > 0, "At least one sample column must be provided"

        conditions = []
        color_position = 0
        sample_map = {}
        color_map = {}

        for sample in sample_cols:
            name_array = sample.split(".")
            replicate = name_array[-1]
            condition = ".".join(name_array[:-1])

            # Add condition if not already present
            if condition not in conditions:
                conditions.append(condition)
                if color_position >= len(payload["settings"]["defaultColorList"]):
                    color_position = 0
                color_map[condition] = payload["settings"]["defaultColorList"][
                    color_position
                ]
                color_position += 1

            # Initialize condition in sample order if needed
            if condition not in payload["settings"]["sampleOrder"]:
                payload["settings"]["sampleOrder"][condition] = []

            # Add sample to its condition
            if sample not in payload["settings"]["sampleOrder"][condition]:
                payload["settings"]["sampleOrder"][condition].append(sample)

            # Make sample visible by default
            if sample not in payload["settings"]["sampleVisible"]:
                payload["settings"]["sampleVisible"][sample] = True

            # Add to sample map
            sample_map[sample] = {
                "condition": condition,
                "replicate": replicate,
                "name": sample,
            }

        payload["settings"]["sampleMap"] = sample_map
        payload["settings"]["colorMap"] = color_map
        payload["settings"]["conditionOrder"] = conditions

    def create_curtain_session_payload(
        self,
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
        primary_id_raw_col: str,
        sample_cols: List[str],
        description: str = "",
    ) -> Dict:
        """
        Create payload for standard Curtain session.

        Args:
            de_file: Path to differential expression file
            raw_file: Path to raw data file
            fc_col: Fold change column name
            transform_fc: Whether to transform fold change values
            transform_significant: Whether to transform significance values
            reverse_fc: Whether to reverse fold change direction
            p_col: Significance column name
            comp_col: Comparison column name (or empty for default)
            comp_select: Comparison selection values
            primary_id_de_col: Primary ID column in differential file
            primary_id_raw_col: Primary ID column in raw file
            sample_cols: Sample column names,
            description: Description

        Returns:
            Configured payload for Curtain session
        """
        return self._prepare_common_payload(
            de_file,
            raw_file,
            fc_col,
            transform_fc,
            transform_significant,
            reverse_fc,
            p_col,
            comp_col,
            comp_select,
            primary_id_de_col,
            primary_id_raw_col,
            sample_cols,
            description,
        )

    def create_curtain_ptm_session_payload(
        self,
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
        primary_id_raw_col: str,
        sample_cols: List[str],
        peptide_col: str,
        acc_col: str,
        position_col: str,
        position_in_peptide_col: str,
        sequence_window_col: str,
        score_col: str,
        description: str = "",
    ) -> Dict:
        """
        Create payload for CurtainPTM session.

        Args:
            de_file: Path to differential expression file
            raw_file: Path to raw data file
            fc_col: Fold change column name
            transform_fc: Whether to transform fold change values
            transform_significant: Whether to transform significance values
            reverse_fc: Whether to reverse fold change direction
            p_col: Significance column name
            comp_col: Comparison column name (or empty for default)
            comp_select: Comparison selection values
            primary_id_de_col: Primary ID column in differential file
            primary_id_raw_col: Primary ID column in raw file
            sample_cols: Sample column names
            peptide_col: Peptide sequence column name
            acc_col: Protein accession column name
            position_col: Position in protein column name
            position_in_peptide_col: Position in peptide column name
            sequence_window_col: Sequence window column name
            score_col: Score column name
            description: Description

        Returns:
            Configured payload for CurtainPTM session
        """
        payload = self._prepare_common_payload(
            de_file,
            raw_file,
            fc_col,
            transform_fc,
            transform_significant,
            reverse_fc,
            p_col,
            comp_col,
            comp_select,
            primary_id_de_col,
            primary_id_raw_col,
            sample_cols,
            description,
        )

        # Add PTM-specific fields
        payload["differentialForm"]["_peptideSequence"] = peptide_col
        payload["differentialForm"]["_accession"] = acc_col
        payload["differentialForm"]["_position"] = position_col
        payload["differentialForm"]["_positionPeptide"] = position_in_peptide_col
        payload["differentialForm"]["_score"] = score_col
        payload["differentialForm"]["_sequenceWindow"] = sequence_window_col

        return payload

    def create_curtain_payload(
        self,
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
        primary_id_raw_col: str,
        sample_cols: List[str],
        peptide_col: str = "",
        acc_col: str = "",
        position_col: str = "",
        position_in_peptide_col: str = "",
        sequence_window_col: str = "",
        score_col: str = "",
    ) -> Dict:
        """
        Create payload for either standard Curtain or CurtainPTM session based on parameters.

        Args:
            de_file: Path to differential expression file
            raw_file: Path to raw data file
            fc_col: Fold change column name
            transform_fc: Whether to transform fold change values
            transform_significant: Whether to transform significance values
            reverse_fc: Whether to reverse fold change direction
            p_col: Significance column name
            comp_col: Comparison column name (or empty for default)
            comp_select: Comparison selection values
            primary_id_de_col: Primary ID column in differential file
            primary_id_raw_col: Primary ID column in raw file
            sample_cols: Sample column names
            peptide_col: Peptide sequence column name (for PTM)
            acc_col: Protein accession column name (for PTM)
            position_col: Position in protein column name (for PTM)
            position_in_peptide_col: Position in peptide column name (for PTM)
            sequence_window_col: Sequence window column name (for PTM)
            score_col: Score column name (for PTM)

        Returns:
            Configured payload for selected Curtain session type
        """
        is_standard = all(
            not col
            for col in [
                peptide_col,
                acc_col,
                position_col,
                position_in_peptide_col,
                sequence_window_col,
                score_col,
            ]
        )

        if is_standard:
            return self.create_curtain_session_payload(
                de_file,
                raw_file,
                fc_col,
                transform_fc,
                transform_significant,
                reverse_fc,
                p_col,
                comp_col,
                comp_select,
                primary_id_de_col,
                primary_id_raw_col,
                sample_cols,
            )
        else:
            return self.create_curtain_ptm_session_payload(
                de_file,
                raw_file,
                fc_col,
                transform_fc,
                transform_significant,
                reverse_fc,
                p_col,
                comp_col,
                comp_select,
                primary_id_de_col,
                primary_id_raw_col,
                sample_cols,
                peptide_col,
                acc_col,
                position_col,
                position_in_peptide_col,
                sequence_window_col,
                score_col,
            )

    def retrieve_curtain_session(self, link_id: str, token: str = "") -> Dict:
        """
        Retrieve Curtain session data.

        Args:
            link_id: Session ID
            token: Optional access token

        Returns:
            Session data

        Raises:
            requests.HTTPError: If API request fails
        """
        url = f"{self.base_url}/curtain/{link_id}/download/token={token}/"
        req = self.request_session.get(url)
        req.raise_for_status()

        res = req.json()
        if "url" in res:
            data_req = self.request_session.get(res["url"])
            data_req.raise_for_status()
            return data_req.json()
        else:
            return res


def parse_curtain_from_json(json_path: str) -> Dict:
    """
    Parse Curtain session data from JSON file.

    Args:
        json_path: Path to JSON file

    Returns:
        Parsed session data
    """
    with open(json_path, "rt") as f:
        data = json.load(f)

    if isinstance(data["settings"], str):
        data["settings"] = json.loads(data["settings"])

    if "version" in data["settings"] and data["settings"]["version"] == "2":
        return parse_v2(data)
    else:
        return parse_old_version(data)


def parse_old_version(json_data: Dict) -> Dict:
    """
    Parse legacy version of Curtain data.

    Args:
        json_data: JSON data to parse

    Returns:
        Processed data
    """
    # Set defaults for missing fields
    if "colormap" not in json_data["settings"]:
        json_data["settings"]["colormap"] = {}
    if "pCutoff" not in json_data["settings"]:
        json_data["settings"]["pCutoff"] = 0.05
    if "log2FCCutoff" not in json_data["settings"]:
        json_data["settings"]["log2FCCutoff"] = 0.6
    if "dataColumns" in json_data["settings"]:
        json_data["settings"]["dataColumns"] = json_data["settings"][
            "dataColumns"
        ].split(",")

    return json_data


def parse_v2(json_data: Dict) -> Dict:
    """
    Parse version 2 format of Curtain data.

    Args:
        json_data: JSON data to parse

    Returns:
        Processed data
    """
    # Process V2 format
    return json_data


def create_imputation_map(
    raw_df: pd.DataFrame | str, primary_id: str, sample_cols: List[str]
) -> Dict:
    if isinstance(raw_df, str):
        if raw_df.endswith("tsv") or raw_df.endswith("txt"):
            raw_df = pd.read_csv(raw_df, sep="\t")
        else:
            raw_df = pd.read_csv(raw_df)
    for sample in sample_cols:
        if sample not in raw_df.columns:
            raise ValueError(f"Sample column '{sample}' not found in raw data")
    if primary_id not in raw_df.columns:
        raise ValueError(f"Primary ID column '{primary_id}' not found in raw data")
    imputation_map = {}
    for i, row in raw_df.iterrows():
        primary_id_value = row[primary_id]
        if primary_id_value not in imputation_map:
            imputation_map[primary_id_value] = {}
        for sample in sample_cols:
            if pd.isnull(row[sample]):
                imputation_map[primary_id_value][sample] = True
    return imputation_map


def add_imputation_map(payload: Dict, imputation_map: Dict) -> Dict:
    """
    Add imputation map to the payload.

    Args:
        payload: Session payload to modify
        imputation_map: Dictionary mapping primary IDs to imputed columns
            The structure should be {primary_id: {column_name: True}}

    Returns:
        Updated payload with imputation map
    """
    payload["settings"]["imputationMap"] = imputation_map
    return payload


def add_uniprot_data(payload: Dict, raw_df: pd.DataFrame | str):
    """
    Add UniProt data to the payload.

    Args:
        payload: Session payload to modify

    Returns:
        Updated payload with UniProt data
    """
    if isinstance(raw_df, str):
        if raw_df.endswith("tsv") or raw_df.endswith("txt"):
            raw_df = pd.read_csv(raw_df, sep="\t")
        else:
            raw_df = pd.read_csv(raw_df)
    parser_columns = "accession,id,gene_names,protein_name,organism_name,organism_id,length,xref_refseq,cc_subcellular_location,sequence,ft_var_seq,cc_alternative_products,cc_function,ft_domain,xref_string,cc_disease,cc_pharmaceutical,ft_mutagen"
    primary_id = payload["rawForm"]["_primaryIDs"]
    if primary_id not in raw_df.columns:
        raise ValueError(f"Primary ID column '{primary_id}' not found in raw data")
    results = {}
    dataMap = {}
    db = {}
    organism = ""
    accMap = {}
    geneNameToAcc = {}
    genesMap = {}
    acc_list = []
    primary_id_map = {}
    primary_id_col = raw_df[primary_id].unique()
    sub_cellular_regex = re.compile(r"[.;]")
    sub_separator_regex = re.compile(r"\s*\{.*?\}\s*")
    domain_position_regex = re.compile(r"(\d+)")
    domain_name_regex = re.compile(r"(.+)")
    uniprot_dataMap = {}
    for pid in primary_id_col:
        if pid not in primary_id_map:
            primary_id_map[pid] = {}
            primary_id_map[pid][pid] = True
        for n in pid.split(";"):
            if pid not in primary_id_map:
                primary_id_map[n] = {}
            primary_id_map[n][pid] = True
        a = pid.split(";")
        dataMap[pid] = pid
        us = UniprotSequence(a[0], True)
        if us.accession:
            if us.accession not in accMap:
                accMap[us.accession] = [us.accession]
            else:
                if us.accession not in accMap[us.accession]:
                    accMap[us.accession].append(us.accession)
            if not us.accession in uniprot_dataMap:
                acc_list.append(us.accession)

    if acc_list:
        parser = UniprotParser(columns=parser_columns)
        db = {}
        allGenes = []
        for res in parser.parse(acc_list, 5000):
            if res:
                df_acc_parse = pd.read_csv(io.StringIO(res), sep="\t")
                organism = df_acc_parse["Organism (ID)"].values[0]
                for i, r in df_acc_parse.iterrows():
                    if pd.notnull(r["Gene Names"]):
                        r["Gene Names"] = r["Gene Names"].replace(" ", ";").upper()
                        df_acc_parse.at[i, "Gene Names"] = r["Gene Names"]
                    if pd.notnull(r["Subcellular location [CC]"]):
                        try:
                            note_position = r["Subcellular location [CC]"].index(
                                "Note="
                            )
                        except ValueError:
                            note_position = -1
                        sub_loc = []
                        if note_position > -1:
                            sub_row = r["Subcellular location [CC]"][:note_position]
                            match = sub_cellular_regex.split(sub_row)
                            if match:
                                for m in match:
                                    if m:
                                        sub_sep_find = sub_separator_regex.sub("", m)
                                        su = sub_sep_find.split(": ")
                                        sub_res = su[-1].strip()
                                        if sub_res:
                                            sub_loc.append(sub_res)
                        df_acc_parse.at[i, "Subcellular location [CC]"] = sub_loc
                    if pd.notnull(r["Domain [FT]"]):
                        domains = []
                        l = 0
                        for d in r["Domain [FT]"].split(";"):
                            if d:
                                if "DOMAIN" in d:
                                    domains.append(dict())
                                    l = len(domains)
                                    match = domain_position_regex.findall(d)
                                    if match:
                                        for m in match:
                                            if m:
                                                if "start" not in domains[l - 1]:
                                                    domains[l - 1]["start"] = int(m)
                                                else:
                                                    domains[l - 1]["end"] = int(m)
                                elif "/note=" in d:
                                    match = domain_name_regex.search(d)
                                    if match:
                                        domains[l - 1]["name"] = match.group(1)
                        df_acc_parse.at[i, "Domain [FT]"] = domains
                    if pd.notnull(r["Mutagenesis"]):
                        mutagenesis = []
                        position = ""
                        for s in r["Mutagenesis"].split("; "):
                            if s:
                                if "MUTAGEN" in s:
                                    position = s.split(" ")[1]
                                elif "/note=" in s:
                                    match = domain_name_regex.search(s)
                                    if match:
                                        mutagenesis.append(
                                            {
                                                "position": position,
                                                "note": match.group(1),
                                            }
                                        )
                        df_acc_parse.at[i, "Mutagenesis"] = mutagenesis
                    df_acc_parse.at[i, "_id"] = r["From"]
                    db[r["Entry"]] = df_acc_parse.iloc[i].fillna("").to_dict()
                    uniprot_dataMap[r["Entry"]] = r["Entry"]
                    uniprot_dataMap[r["From"]] = r["Entry"]
                    if r["Entry"] in accMap:
                        d = accMap[r["Entry"]]
                        for j in d:
                            query = j.replace(",", ";")
                            for q in query.split(";"):
                                if q not in uniprot_dataMap:
                                    uniprot_dataMap[q] = r["Entry"]
                                    if (
                                        pd.notnull(r["Gene Names"])
                                        and r["Gene Names"] != ""
                                    ):
                                        if r["Gene Names"] not in geneNameToAcc:
                                            geneNameToAcc[r["Gene Names"]] = {}
                                        geneNameToAcc[r["Gene Names"]][q] = True
                                else:
                                    if q == r["Entry"]:
                                        uniprot_dataMap[q] = r["Entry"]
                                        if (
                                            pd.notnull(r["Gene Names"])
                                            and r["Gene Names"] != ""
                                        ):
                                            if r["Gene Names"] not in geneNameToAcc:
                                                geneNameToAcc[r["Gene Names"]] = {}
                                            geneNameToAcc[r["Gene Names"]][q] = True
                                    else:
                                        q_splitted = q.split("-")
                                        if len(q_splitted) > 1:
                                            if q_splitted[0] == r["Entry"]:
                                                uniprot_dataMap[q] = r["Entry"]
                                                if (
                                                    pd.notnull(r["Gene Names"])
                                                    and r["Gene Names"] != ""
                                                ):
                                                    if (
                                                        r["Gene Names"]
                                                        not in geneNameToAcc
                                                    ):
                                                        geneNameToAcc[
                                                            r["Gene Names"]
                                                        ] = {}
                                                    geneNameToAcc[r["Gene Names"]][
                                                        q_splitted[0]
                                                    ] = True
    db_df = pd.DataFrame.from_dict(db, orient="index")
    for pid in primary_id_col:
        uniprot = CurtainUniprotData.get_uniprot_data_from_pi_sta(
            pid, accMap, uniprot_dataMap, db_df
        )
        if isinstance(uniprot, pd.Series):
            if pd.notnull(uniprot["Gene Names"]) and uniprot["Gene Names"] != "":
                if uniprot["Gene Names"] not in allGenes:
                    allGenes.append(uniprot["Gene Names"])
                    if uniprot["Gene Names"] not in genesMap:
                        genesMap[uniprot["Gene Names"]] = {}
                        genesMap[uniprot["Gene Names"]][uniprot["Gene Names"]] = True
                    for n in uniprot["Gene Names"].split(";"):
                        if n not in genesMap:
                            genesMap[n] = {}
                        genesMap[n][uniprot["Gene Names"]] = True

    if "extraData" not in payload:
        payload["extraData"] = {
            "uniprot": {
                "accMap": {"dataType": "Map", "value": list(accMap.items())},
                "dataMap": {"dataType": "Map", "value": list(uniprot_dataMap.items())},
                "db": {"dataType": "Map", "value": list(db.items())},
                "organism": organism,
                "geneNameToAcc": geneNameToAcc,
                "results": {"dataType": "Map", "value": list(results.items())},
            },
            "data": {
                "dataMap": {"dataType": "Map", "value": list(dataMap.items())},
                "genesMap": genesMap,
                "primaryIDsmap": primary_id_map,
                "allGenes": allGenes,
            },
        }
    
    # Set fetchUniProt to indicate UniProt data has been added
    payload["fetchUniProt"] = True


def add_uniprot_data_ptm(payload: Dict, raw_df: pd.DataFrame | str, de_df: pd.DataFrame | str):
    """
    Add UniProt data to CurtainPTM payload with PTM-specific mappings.

    Args:
        payload: CurtainPTM session payload to modify
        raw_df: Raw data DataFrame or file path
        de_df: Differential expression DataFrame or file path (needed for accession mapping)

    Returns:
        Updated payload with CurtainPTM-specific UniProt data and mappings
    """
    # Load DataFrames if paths provided
    if isinstance(raw_df, str):
        if raw_df.endswith("tsv") or raw_df.endswith("txt"):
            raw_df = pd.read_csv(raw_df, sep="\t")
        else:
            raw_df = pd.read_csv(raw_df)
    
    if isinstance(de_df, str):
        if de_df.endswith("tsv") or de_df.endswith("txt"):
            de_df = pd.read_csv(de_df, sep="\t")
        else:
            de_df = pd.read_csv(de_df)
    
    # Enhanced UniProt fields for PTM analysis
    parser_columns = "accession,id,gene_names,protein_name,organism_name,organism_id,length,cc_subcellular_location,sequence,ft_var_seq,cc_alternative_products,ft_domain,xref_string,ft_mod_res,cc_function,cc_disease,cc_pharmaceutical,ft_mutagen,xref_mim"
    
    # Get column mappings
    primary_id_raw = payload["rawForm"]["_primaryIDs"]
    primary_id_de = payload["differentialForm"]["_primaryIDs"]
    accession_col = payload["differentialForm"]["_accession"]
    
    # Validate columns exist
    if primary_id_raw not in raw_df.columns:
        raise ValueError(f"Primary ID column '{primary_id_raw}' not found in raw data")
    if primary_id_de not in de_df.columns:
        raise ValueError(f"Primary ID column '{primary_id_de}' not found in differential data")
    if accession_col not in de_df.columns:
        raise ValueError(f"Accession column '{accession_col}' not found in differential data")
    
    # Initialize mapping structures
    results = {}
    dataMap = {}
    db = {}
    organism = ""
    accMap = {}
    geneNameToPrimary = {}  # CurtainPTM uses geneNameToPrimary instead of geneNameToAcc
    genesMap = {}
    
    # CurtainPTM-specific structures
    accessionMap = {}
    accessionToPrimaryIDs = {}
    accessionList = []
    primaryIDsList = []
    
    acc_list = []
    uniprot_dataMap = {}
    
    # Build primary ID mappings from raw data
    primary_id_col_raw = raw_df[primary_id_raw].unique()
    primaryIDsList = list(primary_id_col_raw)
    
    # Build accession mappings from differential data
    accession_col_data = de_df[accession_col].dropna().unique()
    accessionList = list(accession_col_data)
    
    # Build accession map (PTM-specific)
    for accession in accessionList:
        if accession not in accessionMap:
            accessionMap[accession] = {}
            accessionMap[accession][accession] = True
        
        # Handle semicolon-separated accessions
        for split_acc in accession.split(";"):
            if split_acc not in accessionMap:
                accessionMap[split_acc] = {}
            accessionMap[split_acc][accession] = True
    
    # Build accessionToPrimaryIDs mapping
    for i, row in de_df.iterrows():
        accession = str(row[accession_col])
        primary_id = str(row[primary_id_de])
        
        if pd.notna(accession) and pd.notna(primary_id):
            accession_parts = accession.split(";")
            us = UniprotSequence(accession_parts[0], True)
            
            if us.accession:
                if us.accession not in accessionToPrimaryIDs:
                    accessionToPrimaryIDs[us.accession] = {}
                accessionToPrimaryIDs[us.accession][primary_id] = True
                
                # Build dataMap entries
                dataMap[accession] = accession
                dataMap[primary_id] = accession
                
                # Prepare for UniProt query
                if us.accession not in accMap:
                    accMap[us.accession] = [us.accession]
                else:
                    if us.accession not in accMap[us.accession]:
                        accMap[us.accession].append(us.accession)
                
                if us.accession not in uniprot_dataMap:
                    acc_list.append(us.accession)
    
    # Process UniProt data if accessions found
    if acc_list:
        parser = UniprotParser(columns=parser_columns)
        allGenes = []
        
        for res in parser.parse(acc_list, 5000):
            if res:
                df_acc_parse = pd.read_csv(io.StringIO(res), sep="\t")
                organism = df_acc_parse["Organism (ID)"].values[0]
                
                for i, r in df_acc_parse.iterrows():
                    # Process gene names
                    if pd.notnull(r["Gene Names"]):
                        r["Gene Names"] = r["Gene Names"].replace(" ", ";").upper()
                        df_acc_parse.at[i, "Gene Names"] = r["Gene Names"]
                    
                    # Process subcellular location
                    if pd.notnull(r["Subcellular location [CC]"]):
                        try:
                            note_position = r["Subcellular location [CC]"].index("Note=")
                        except ValueError:
                            note_position = -1
                        
                        sub_loc = []
                        if note_position > -1:
                            sub_row = r["Subcellular location [CC]"][:note_position]
                            match = re.split(r"[.;]", sub_row)
                            if match:
                                for m in match:
                                    if m:
                                        sub_sep_find = re.sub(r"\s*\{.*?\}\s*", "", m)
                                        su = sub_sep_find.split(": ")
                                        sub_res = su[-1].strip()
                                        if sub_res:
                                            sub_loc.append(sub_res)
                        df_acc_parse.at[i, "Subcellular location [CC]"] = sub_loc
                    
                    # Process modified residues (PTM-specific)
                    if pd.notnull(r["Modified residue"]):
                        mod_res = []
                        mods = r["Modified residue"].split("; ")
                        mod_position = -1
                        mod_type = ""
                        
                        for m in mods:
                            if m.startswith("MOD_RES"):
                                if ":" in m:
                                    mod_position = int(m.split(":")[1]) - 1
                                else:
                                    mod_position = int(m.split(" ")[1]) - 1
                            elif "note=" in m:
                                modre = re.search(r'".+"', m)
                                if modre:
                                    mod_type = modre.group(0)
                                    if mod_position >= 0 and mod_position < len(r["Sequence"]):
                                        mod_res.append({
                                            "position": mod_position + 1,
                                            "residue": r["Sequence"][mod_position],
                                            "modType": mod_type.replace('"', "")
                                        })
                        df_acc_parse.at[i, "Modified residue"] = mod_res
                    
                    # Process domains
                    if pd.notnull(r["Domain [FT]"]):
                        domains = []
                        l = 0
                        for d in r["Domain [FT]"].split(";"):
                            if d:
                                if "DOMAIN" in d:
                                    domains.append({})
                                    l = len(domains)
                                    match = re.findall(r"(\d+)", d)
                                    if match:
                                        for m in match:
                                            if m:
                                                if "start" not in domains[l - 1]:
                                                    domains[l - 1]["start"] = int(m)
                                                else:
                                                    domains[l - 1]["end"] = int(m)
                                elif "/note=" in d:
                                    match = re.search(r'"(.+)"', d)
                                    if match:
                                        domains[l - 1]["name"] = match.group(1)
                        df_acc_parse.at[i, "Domain [FT]"] = domains
                    
                    # Process mutagenesis
                    if pd.notnull(r["Mutagenesis"]):
                        mutagenesis = []
                        position = ""
                        for s in r["Mutagenesis"].split("; "):
                            if s:
                                if "MUTAGEN" in s:
                                    position = s.split(" ")[1]
                                elif "/note=" in s:
                                    match = re.search(r'"(.+)"', s)
                                    if match:
                                        mutagenesis.append({
                                            "position": position,
                                            "note": match.group(1)
                                        })
                        df_acc_parse.at[i, "Mutagenesis"] = mutagenesis
                    
                    df_acc_parse.at[i, "_id"] = r["From"]
                    db[r["Entry"]] = df_acc_parse.iloc[i].fillna("").to_dict()
                    uniprot_dataMap[r["Entry"]] = r["Entry"]
                    uniprot_dataMap[r["From"]] = r["Entry"]
                    
                    # Build gene name to primary ID mapping (CurtainPTM-specific)
                    if pd.notnull(r["Gene Names"]) and r["Gene Names"] != "":
                        gene_name = r["Gene Names"]
                        if gene_name not in allGenes:
                            allGenes.append(gene_name)
                        
                        if gene_name not in genesMap:
                            genesMap[gene_name] = {}
                            genesMap[gene_name][gene_name] = True
                        
                        for n in gene_name.split(";"):
                            if n not in genesMap:
                                genesMap[n] = {}
                            genesMap[n][gene_name] = True
                        
                        # Map gene names to primary IDs
                        if gene_name not in geneNameToPrimary:
                            geneNameToPrimary[gene_name] = {}
                        
                        if r["Entry"] in accessionToPrimaryIDs:
                            for primary_id in accessionToPrimaryIDs[r["Entry"]]:
                                geneNameToPrimary[gene_name][primary_id] = True
                    
                    # Handle accession mappings
                    if r["Entry"] in accMap:
                        d = accMap[r["Entry"]]
                        for j in d:
                            query = j.replace(",", ";")
                            for q in query.split(";"):
                                if q not in uniprot_dataMap:
                                    uniprot_dataMap[q] = r["Entry"]
                                else:
                                    if q == r["Entry"]:
                                        uniprot_dataMap[q] = r["Entry"]
                                    else:
                                        q_splitted = q.split("-")
                                        if len(q_splitted) > 1:
                                            if q_splitted[0] == r["Entry"]:
                                                uniprot_dataMap[q] = r["Entry"]
    
    # Build the CurtainPTM-specific extraData structure
    payload["extraData"] = {
        "uniprot": {
            "accMap": {"dataType": "Map", "value": list(accMap.items())},
            "dataMap": {"dataType": "Map", "value": list(uniprot_dataMap.items())},
            "db": {"dataType": "Map", "value": list(db.items())},
            "organism": organism,
            "geneNameToPrimary": geneNameToPrimary,  # CurtainPTM uses this instead of geneNameToAcc
            "results": {"dataType": "Map", "value": list(results.items())},
        },
        "data": {
            # CurtainPTM-specific data structures
            "accessionToPrimaryIDs": accessionToPrimaryIDs,
            "primaryIDsList": primaryIDsList,
            "accessionList": accessionList,
            "accessionMap": accessionMap,
            "genesMap": genesMap,
            "allGenes": allGenes,
            "dataMap": {"dataType": "Map", "value": list(dataMap.items())},
        },
    }
    
    # Set fetchUniProt to indicate UniProt data has been added
    payload["fetchUniProt"] = True
    
    return payload


def configure_volcano_plot(payload: Dict, **kwargs) -> Dict:
    """
    Configure volcano plot settings in payload.
    
    Available parameters:
    - x_min, x_max, y_min, y_max: Axis ranges
    - x_title, y_title: Axis titles (default: "Log2FC", "-log10(p-value)")
    - x_tick_interval, y_tick_interval: Tick intervals (dtickX, dtickY)
    - x_tick_length, y_tick_length: Tick lengths (default: 5)
    - width, height: Plot dimensions (default: 800, 1000)
    - margin_left, margin_right, margin_bottom, margin_top: Plot margins
    - show_x_grid, show_y_grid: Grid line visibility (default: True)
    - title: Plot title
    - legend_x, legend_y: Legend position
    - marker_size: Scatter plot marker size (default: 10)
    """
    if 'settings' not in payload:
        payload['settings'] = {}
    
    # Initialize volcanoAxis if not exists
    if 'volcanoAxis' not in payload['settings']:
        payload['settings']['volcanoAxis'] = {
            'minX': None, 'maxX': None, 'minY': None, 'maxY': None,
            'x': "Log2FC", 'y': "-log10(p-value)",
            'dtickX': None, 'dtickY': None,
            'ticklenX': 5, 'ticklenY': 5
        }
    
    # Axis ranges
    if 'x_min' in kwargs:
        payload['settings']['volcanoAxis']['minX'] = kwargs['x_min']
    if 'x_max' in kwargs:
        payload['settings']['volcanoAxis']['maxX'] = kwargs['x_max']
    if 'y_min' in kwargs:
        payload['settings']['volcanoAxis']['minY'] = kwargs['y_min']
    if 'y_max' in kwargs:
        payload['settings']['volcanoAxis']['maxY'] = kwargs['y_max']
    
    # Axis titles
    if 'x_title' in kwargs:
        payload['settings']['volcanoAxis']['x'] = kwargs['x_title']
    if 'y_title' in kwargs:
        payload['settings']['volcanoAxis']['y'] = kwargs['y_title']
    
    # Tick settings
    if 'x_tick_interval' in kwargs:
        payload['settings']['volcanoAxis']['dtickX'] = kwargs['x_tick_interval']
    if 'y_tick_interval' in kwargs:
        payload['settings']['volcanoAxis']['dtickY'] = kwargs['y_tick_interval']
    if 'x_tick_length' in kwargs:
        payload['settings']['volcanoAxis']['ticklenX'] = kwargs['x_tick_length']
    if 'y_tick_length' in kwargs:
        payload['settings']['volcanoAxis']['ticklenY'] = kwargs['y_tick_length']
    
    # Plot dimensions
    if 'volcanoPlotDimension' not in payload['settings']:
        payload['settings']['volcanoPlotDimension'] = {
            'width': 800, 'height': 1000,
            'margin': {'l': None, 'r': None, 'b': None, 't': None}
        }
    
    if 'width' in kwargs:
        payload['settings']['volcanoPlotDimension']['width'] = kwargs['width']
    if 'height' in kwargs:
        payload['settings']['volcanoPlotDimension']['height'] = kwargs['height']
    
    # Margins
    if any(k in kwargs for k in ['margin_left', 'margin_right', 'margin_bottom', 'margin_top']):
        if 'margin_left' in kwargs:
            payload['settings']['volcanoPlotDimension']['margin']['l'] = kwargs['margin_left']
        if 'margin_right' in kwargs:
            payload['settings']['volcanoPlotDimension']['margin']['r'] = kwargs['margin_right']
        if 'margin_bottom' in kwargs:
            payload['settings']['volcanoPlotDimension']['margin']['b'] = kwargs['margin_bottom']
        if 'margin_top' in kwargs:
            payload['settings']['volcanoPlotDimension']['margin']['t'] = kwargs['margin_top']
    
    # Grid lines
    if 'volcanoPlotGrid' not in payload['settings']:
        payload['settings']['volcanoPlotGrid'] = {'x': True, 'y': True}
    
    if 'show_x_grid' in kwargs:
        payload['settings']['volcanoPlotGrid']['x'] = kwargs['show_x_grid']
    if 'show_y_grid' in kwargs:
        payload['settings']['volcanoPlotGrid']['y'] = kwargs['show_y_grid']
    
    # Title and legend
    if 'title' in kwargs:
        payload['settings']['volcanoPlotTitle'] = kwargs['title']
    if 'legend_x' in kwargs:
        payload['settings']['volcanoPlotLegendX'] = kwargs['legend_x']
    if 'legend_y' in kwargs:
        payload['settings']['volcanoPlotLegendY'] = kwargs['legend_y']
    
    # Marker size
    if 'marker_size' in kwargs:
        payload['settings']['scatterPlotMarkerSize'] = kwargs['marker_size']
    
    # Additional shapes (for custom annotations)
    if 'additional_shapes' in kwargs:
        payload['settings']['volcanoAdditionalShapes'] = kwargs['additional_shapes']
    
    return payload


def configure_bar_chart(payload: Dict, **kwargs) -> Dict:
    """
    Configure bar chart settings in payload.
    
    Available parameters:
    - bar_chart_width: Width per column in individual bar chart
    - average_bar_chart_width: Width per column in average bar chart
    - violin_plot_width: Width per column in violin plot
    - profile_plot_width: Width per column in profile plot
    - condition_colors: Dictionary mapping condition names to colors (overrides default colorMap)
    - violin_point_position: Position of violin plot points (-2 = hide, default: -2)
    """
    if 'settings' not in payload:
        payload['settings'] = {}
    
    # Initialize columnSize if not exists
    if 'columnSize' not in payload['settings']:
        payload['settings']['columnSize'] = {
            'barChart': 0, 'averageBarChart': 0,
            'violinPlot': 0, 'profilePlot': 0
        }
    
    # Column sizing
    if 'bar_chart_width' in kwargs:
        payload['settings']['columnSize']['barChart'] = kwargs['bar_chart_width']
    if 'average_bar_chart_width' in kwargs:
        payload['settings']['columnSize']['averageBarChart'] = kwargs['average_bar_chart_width']
    if 'violin_plot_width' in kwargs:
        payload['settings']['columnSize']['violinPlot'] = kwargs['violin_plot_width']
    if 'profile_plot_width' in kwargs:
        payload['settings']['columnSize']['profilePlot'] = kwargs['profile_plot_width']
    
    # Bar chart specific color overrides
    if 'condition_colors' in kwargs:
        if 'barchartColorMap' not in payload['settings']:
            payload['settings']['barchartColorMap'] = {}
        payload['settings']['barchartColorMap'].update(kwargs['condition_colors'])
    
    # Violin plot settings
    if 'violin_point_position' in kwargs:
        payload['settings']['violinPointPos'] = kwargs['violin_point_position']
    
    return payload


def configure_general_plot_settings(payload: Dict, **kwargs) -> Dict:
    """
    Configure general plot settings that affect all visualizations.
    
    Available parameters:
    - font_family: Font family for all plots (default: "Arial")
    - p_cutoff: P-value significance cutoff (default: 0.05)
    - fc_cutoff: Log2 fold change cutoff (default: 0.6)
    - default_colors: List of default colors for conditions
    - condition_colors: Dictionary mapping condition names to colors
    - condition_order: List specifying order of conditions in plots
    - sample_visibility: Dictionary controlling which samples are visible
    """
    if 'settings' not in payload:
        payload['settings'] = {}
    
    # Font settings
    if 'font_family' in kwargs:
        payload['settings']['plotFontFamily'] = kwargs['font_family']
    
    # Significance thresholds
    if 'p_cutoff' in kwargs:
        payload['settings']['pCutoff'] = kwargs['p_cutoff']
    if 'fc_cutoff' in kwargs:
        payload['settings']['log2FCCutoff'] = kwargs['fc_cutoff']
    
    # Color settings
    if 'default_colors' in kwargs:
        payload['settings']['defaultColorList'] = kwargs['default_colors']
    
    if 'condition_colors' in kwargs:
        if 'colorMap' not in payload['settings']:
            payload['settings']['colorMap'] = {}
        payload['settings']['colorMap'].update(kwargs['condition_colors'])
    
    # Sample and condition management
    if 'condition_order' in kwargs:
        payload['settings']['conditionOrder'] = kwargs['condition_order']
    
    if 'sample_visibility' in kwargs:
        if 'sampleVisible' not in payload['settings']:
            payload['settings']['sampleVisible'] = {}
        payload['settings']['sampleVisible'].update(kwargs['sample_visibility'])
    
    return payload


def configure_ptm_specific_settings(payload: Dict, **kwargs) -> Dict:
    """
    Configure PTM-specific settings (CurtainPTM only).
    
    Available parameters:
    - custom_ptm_data: Custom PTM database annotations in the format of a nested dictionary where the first level keys are the name of the database, the second level is the uniprot accession id, and the third level is the primary ID and the last value is a list of dictionaries with {position: int, residue: str} where position is 0-based index.
    - variant_corrections: PTM position variant corrections is a dictionary where the keys are the uniprot accession IDs and the values are the actual uniprot accession IDs of those accession IDs found in the datasets
    - custom_sequences: Custom peptide sequence definitions
    """
    if 'settings' not in payload:
        payload['settings'] = {}
    
    # PTM-specific settings
    if 'custom_ptm_data' in kwargs:
        payload['settings']['customPTMData'] = kwargs['custom_ptm_data']
    
    if 'variant_corrections' in kwargs:
        payload['settings']['variantCorrection'] = kwargs['variant_corrections']
    
    if 'custom_sequences' in kwargs:
        payload['settings']['customSequences'] = kwargs['custom_sequences']
    
    return payload
