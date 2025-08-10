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
    """Custom JSON encoder for handling NumPy data types.
    
    This encoder extends the standard JSON encoder to handle NumPy-specific
    data types that are commonly used in scientific computing.
    """
    
    def default(self, obj):
        """Override default JSON encoder to handle NumPy types.
        
        Args:
            obj: Object to be JSON-serialized
            
        Returns:
            JSON-serializable representation of the object
        """
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
        """Convert database dictionary to DataFrame.
        
        Args:
            db: Dictionary containing database values with 'value' key
            
        Returns:
            DataFrame constructed from the database values
        """
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
        """Get UniProt data for a primary ID using static lookup.
        
        This is a static version of get_uniprot_data_from_pi that can be used
        without instantiating the CurtainUniprotData class.
        
        Args:
            primary_id: Primary identifier to search for
            accMap: Mapping of primary IDs to accession numbers
            dataMap: Mapping of accessions to internal IDs
            db: DataFrame containing UniProt database entries
            
        Returns:
            DataFrame row containing the UniProt data or None if not found
        """
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
        """Get authentication headers if API key is provided.
        
        Returns:
            Dictionary containing authentication headers or empty dict if no API key
        """
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

                                                        if pbar:
                                                            pbar.update(len(chunk))

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
                                            time.sleep(2**dl_attempt)
                                            continue
                                        else:
                                            return None
                                else:
                                    result = requests.get(data["url"])
                                    if result.status_code == 200:
                                        return result.json()
                                    elif dl_attempt < retries:
                                        time.sleep(2**dl_attempt)
                                        continue
                                    else:
                                        return None

                            except (
                                requests.exceptions.RequestException,
                                json.JSONDecodeError,
                            ) as e:
                                if dl_attempt < retries:
                                    time.sleep(2**dl_attempt)
                                    continue
                                return None
                    else:
                        return data
                elif attempt < retries:
                    time.sleep(2**attempt)
                    continue
                else:
                    return None

            except requests.exceptions.RequestException as e:
                if attempt < retries:
                    time.sleep(2**attempt)
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

        with open(de_file, "rt") as f, open(raw_file, "rt") as f2:
            payload["processed"] = f.read()
            payload["raw"] = f2.read()
        payload["settings"]["description"] = description
        
        payload["differentialForm"]["_foldChange"] = fc_col
        payload["differentialForm"]["_significant"] = p_col
        payload["differentialForm"]["_comparison"] = comp_col or "CurtainSetComparison"
        payload["differentialForm"]["_comparisonSelect"] = comp_select or ["1"]
        payload["differentialForm"]["_primaryIDs"] = primary_id_de_col
        payload["rawForm"]["_primaryIDs"] = primary_id_raw_col
        payload["rawForm"]["_samples"] = sample_cols

        payload["differentialForm"]["_transformFC"] = transform_fc
        payload["differentialForm"]["_transformSignificant"] = transform_significant
        payload["differentialForm"]["_reverseFoldChange"] = reverse_fc

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

            if condition not in conditions:
                conditions.append(condition)
                if color_position >= len(payload["settings"]["defaultColorList"]):
                    color_position = 0
                color_map[condition] = payload["settings"]["defaultColorList"][
                    color_position
                ]
                color_position += 1

            if condition not in payload["settings"]["sampleOrder"]:
                payload["settings"]["sampleOrder"][condition] = []

            if sample not in payload["settings"]["sampleOrder"][condition]:
                payload["settings"]["sampleOrder"][condition].append(sample)

            if sample not in payload["settings"]["sampleVisible"]:
                payload["settings"]["sampleVisible"][sample] = True

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


class Payload:
    """
    Object-oriented interface for building and configuring Curtain/CurtainPTM payloads.
    
    This class provides a clean, chainable API for creating and configuring
    data submissions to Curtain and CurtainPTM platforms. It wraps the existing
    functional API while maintaining state consistency.
    
    Example:
        >>> payload = Payload(curtain_type='TP')
        >>> payload.configure_volcano_plot(
        ...     x_title="Log2 Fold Change",
        ...     y_title="-log10(p-value)"
        ... ).add_search_group(
        ...     "Kinases", 
        ...     gene_names=["AKT1", "MTOR"]
        ... )
    """
    
    def __init__(self, curtain_type: str = 'TP', initial_payload: Dict = None):
        """
        Initialize a new Payload instance.
        
        Args:
            curtain_type: Type of Curtain platform ('TP' or 'PTM')
            initial_payload: Existing payload data to start from
        """
        if initial_payload is not None:
            self.data = deepcopy(initial_payload)
        else:
            self.data = deepcopy(curtain_base_payload)
        
        # Validate curtain type
        if curtain_type.upper() not in ['TP', 'PTM']:
            raise ValueError("curtain_type must be 'TP' or 'PTM'")
        
        self.curtain_type = curtain_type.upper()
        self.data['curtain_type'] = self.curtain_type
        
        # Ensure required structures exist
        self._ensure_structure()
    
    def _ensure_structure(self):
        """Ensure all required payload structures exist."""
        required_keys = {
            'settings': {},
            'selections': [],
            'selectionsMap': {},
            'selectionsName': [],
            'extraData': {}
        }
        
        for key, default_value in required_keys.items():
            if key not in self.data:
                self.data[key] = default_value
        
        if 'colorMap' not in self.data['settings']:
            self.data['settings']['colorMap'] = {}
    
    @property
    def is_curtain_ptm(self) -> bool:
        """Check if this payload is for CurtainPTM."""
        return self.curtain_type == 'PTM'
    
    def to_dict(self) -> Dict:
        """Return the payload data as a dictionary."""
        return deepcopy(self.data)
    
    def to_json(self, **kwargs) -> str:
        """Return the payload data as JSON string."""
        return json.dumps(self.data, cls=NumpyEncoder, **kwargs)
    
    # Configuration methods - properly delegating to module-level functions
    def configure_volcano_plot(self, **kwargs) -> 'Payload':
        """Configure volcano plot settings. Returns self for method chaining."""
        self.data = configure_volcano_plot(self.data, **kwargs)
        return self
    
    def configure_bar_chart(self, **kwargs) -> 'Payload':
        """Configure bar chart settings. Returns self for method chaining."""
        self.data = configure_bar_chart(self.data, **kwargs)
        return self
    
    def configure_general_plot_settings(self, **kwargs) -> 'Payload':
        """Configure general plot settings. Returns self for method chaining."""
        self.data = configure_general_plot_settings(self.data, **kwargs)
        return self
    
    def configure_ptm_specific_settings(self, **kwargs) -> 'Payload':
        """Configure PTM-specific settings. Returns self for method chaining."""
        self.data = configure_ptm_specific_settings(self.data, **kwargs)
        return self
    
    def configure_sample_conditions(self, sample_condition_map: Dict[str, str], 
                                   condition_colors: Dict[str, str] = None, 
                                   condition_order: List[str] = None) -> 'Payload':
        """Configure custom sample-to-condition mapping. Returns self for method chaining."""
        self.data = configure_sample_conditions(self.data, sample_condition_map, condition_colors, condition_order)
        return self
    
    def add_search_group(self, group_name: str, protein_ids: List[str] = None,
                        gene_names: List[str] = None, color: str = None, 
                        specific_comparison_label: str = None, raw_df: pd.DataFrame = None,
                        primary_id_col: str = None, gene_name_col: str = None) -> 'Payload':
        """Add a search group. Returns self for method chaining."""
        self.data = add_search_group_with_type(self.data, group_name, self.curtain_type, protein_ids, gene_names, color, 
                                    specific_comparison_label, raw_df, primary_id_col, gene_name_col)
        return self
    
    def create_batch_search_group(self, group_name: str, search_input: str, 
                                 search_type: str = "Gene Names", color: str = None,
                                 specific_comparison_label: str = None) -> 'Payload':
        """Create search group from batch input (like frontend batch search). Returns self for method chaining."""
        self.data = create_batch_search_group(self.data, group_name, search_input, search_type, 
                                            color, specific_comparison_label)
        return self
    
    def add_multiple_search_groups(self, groups: Dict[str, Dict]) -> 'Payload':
        """Add multiple search groups at once. Returns self for method chaining."""
        self.data = add_multiple_search_groups(self.data, groups)
        return self
    
    def create_search_group_from_regex(self, group_name: str, pattern: str,
                                      raw_df: pd.DataFrame, search_column: str,
                                      primary_id_col: str, color: str = None,
                                      specific_comparison_label: str = None) -> 'Payload':
        """Create search group using regex pattern. Returns self for method chaining."""
        self.data = create_search_group_from_regex(self.data, group_name, pattern, raw_df, 
                                                  search_column, primary_id_col, color, specific_comparison_label)
        return self
    
    def remove_search_group(self, group_name: str) -> 'Payload':
        """Remove a search group. Returns self for method chaining."""
        self.data = remove_search_group(self.data, group_name)
        return self
    
    def clear_search_groups(self) -> 'Payload':
        """Clear all search groups. Returns self for method chaining."""
        if 'selections' in self.data:
            self.data['selections'] = []
        if 'selectionsMap' in self.data:
            self.data['selectionsMap'] = {}
        if 'selectionsName' in self.data:
            self.data['selectionsName'] = []
        return self
    
    def get_search_groups(self) -> List[str]:
        """Get list of search group names."""
        return self.data.get('selectionsName', [])
    
    def get_conditions(self) -> List[str]:
        """Get list of condition names from sample mapping."""
        if 'settings' in self.data and 'sampleMap' in self.data['settings']:
            sample_map = self.data['settings']['sampleMap']
            return list(set(sample_map[s]['condition'] for s in sample_map if 'condition' in sample_map[s]))
        return []
    
    def create_significance_groups(self, p_cutoff: float = 0.05, fc_cutoff: float = 0.6, 
                                 colors: Dict[str, str] = None) -> 'Payload':
        """
        Create automatic significance groups based on p-value and fold change cutoffs.
        
        This replicates the frontend's automatic grouping system that categorizes
        datapoints into 4 groups based on significance and fold change thresholds.
        
        Args:
            p_cutoff: P-value threshold for significance (default: 0.05).
            fc_cutoff: Log2 fold change threshold (default: 0.6).
            colors: Optional dictionary mapping group types to colors.
                   Keys should be like "P-value > ;FC > " format.
        
        Returns:
            Self for method chaining.
        """
        if 'settings' not in self.data:
            self.data['settings'] = {}
        
        self.data['settings']['pCutoff'] = p_cutoff
        self.data['settings']['log2FCCutoff'] = fc_cutoff
        
        if colors:
            if 'colorMap' not in self.data['settings']:
                self.data['settings']['colorMap'] = {}
            
            self._apply_significance_colors(colors, p_cutoff, fc_cutoff)
        
        return self
    
    def get_significance_group(self, x: float, y: float, comparison_label: str = None) -> tuple:
        """
        Determine significance group for a datapoint based on x,y coordinates.
        
        Replicates the frontend significantGroup() function logic with comparison labels.
        
        Args:
            x: Log2 fold change value.
            y: -log10(p-value) value.
            comparison_label: Comparison label to append to group name. If None, uses first available comparison.
        
        Returns:
            Tuple of (full_group_name, position_key) where:
            - full_group_name: Complete group name with comparison label suffix
            - position_key: Position identifier for consistent coloring
        """
        p_cutoff = self.data.get('settings', {}).get('pCutoff', 0.05)
        fc_cutoff = self.data.get('settings', {}).get('log2FCCutoff', 0.6)
        
        import math
        ylog = -math.log10(p_cutoff)
        groups = []
        position = ""
        
        if ylog > y:
            groups.append(f"P-value > {p_cutoff}")
            position = "P-value > "
        else:
            groups.append(f"P-value <= {p_cutoff}")
            position = "P-value <= "
        
        if abs(x) > fc_cutoff:
            groups.append(f"FC > {fc_cutoff}")
            position += "FC > "
        else:
            groups.append(f"FC <= {fc_cutoff}")
            position += "FC <= "
        
        # Get comparison label for group naming
        if comparison_label is None:
            diff_form = self.data.get('differentialForm', {})
            comparison_select = diff_form.get('_comparisonSelect', ['1'])
            if isinstance(comparison_select, list) and comparison_select:
                comparison_label = comparison_select[0]
            else:
                comparison_label = '1'
        
        base_group_name = ";".join(groups)
        
        # Add comparison suffix based on curtain type
        if self.curtain_type == 'PTM':
            # PTM groups don't have comparison suffix
            full_group_name = base_group_name
        else:
            # TP groups have comparison suffix
            full_group_name = f"{base_group_name} ({comparison_label})"
        
        return (full_group_name, position)
    
    def _apply_significance_colors(self, colors: Dict[str, str], p_cutoff: float, fc_cutoff: float):
        """Apply colors to significance group patterns."""
        diff_form = self.data.get('differentialForm', {})
        comparison_select = diff_form.get('_comparisonSelect', ['1'])
        
        if isinstance(comparison_select, str):
            comparison_select = [comparison_select]
        
        for comparison in comparison_select:
            for position_key, color in colors.items():
                # Build the actual group name with values
                parts = []
                position_parts = position_key.split('FC')
                
                # Handle P-value part
                if 'P-value > ' in position_parts[0]:
                    parts.append(f"P-value > {p_cutoff}")
                elif 'P-value <= ' in position_parts[0]:
                    parts.append(f"P-value <= {p_cutoff}")
                
                # Handle FC part  
                if len(position_parts) > 1:
                    fc_part = position_parts[1].strip()
                    if fc_part.startswith('> '):
                        parts.append(f"FC > {fc_cutoff}")
                    elif fc_part.startswith('<= '):
                        parts.append(f"FC <= {fc_cutoff}")
                
                base_group_name = ";".join(parts)
                
                # Add comparison suffix based on curtain type
                if self.curtain_type == 'PTM':
                    full_group_name = base_group_name
                else:
                    full_group_name = f"{base_group_name} ({comparison})"
                
                self.data['settings']['colorMap'][full_group_name] = color
    
    def set_significance_group_colors(self, colors: Dict[str, str]) -> 'Payload':
        """
        Set colors for significance groups.
        
        Args:
            colors: Dictionary mapping significance group patterns to hex colors.
                   Keys should be position patterns like:
                   - "P-value > FC > " (non-significant, high fold change)
                   - "P-value > FC <= " (non-significant, low fold change)  
                   - "P-value <= FC > " (significant, high fold change)
                   - "P-value <= FC <= " (significant, low fold change)
        
        Returns:
            Self for method chaining.
        """
        if 'settings' not in self.data:
            self.data['settings'] = {}
        if 'colorMap' not in self.data['settings']:
            self.data['settings']['colorMap'] = {}
        
        p_cutoff = self.data.get('settings', {}).get('pCutoff', 0.05)
        fc_cutoff = self.data.get('settings', {}).get('log2FCCutoff', 0.6)
        
        self._apply_significance_colors(colors, p_cutoff, fc_cutoff)
        
        return self
    
    def get_default_significance_colors(self) -> Dict[str, str]:
        """
        Get the default color mapping for significance groups.
        
        Returns:
            Dictionary mapping significance group patterns to default colors.
        """
        default_colors = [
            "#fd7f6f",  # Light red
            "#7eb0d5",  # Light blue  
            "#b2e061",  # Light green
            "#bd7ebe",  # Light purple
        ]
        
        patterns = [
            "P-value > FC <= ",   # Non-significant, low FC - typically grey/muted
            "P-value > FC > ",    # Non-significant, high FC - typically orange
            "P-value <= FC <= ",  # Significant, low FC - typically blue
            "P-value <= FC > ",   # Significant, high FC - typically red (most interesting)
        ]
        
        return dict(zip(patterns, default_colors))
    
    def get_samples(self) -> List[str]:
        """Get list of sample names."""
        if 'settings' in self.data and 'sampleMap' in self.data['settings']:
            return list(self.data['settings']['sampleMap'].keys())
        return []
    
    def detect_sample_patterns(self, raw_df: pd.DataFrame, sample_columns: List[str]) -> Dict[str, str]:
        """Detect sample condition patterns from column names."""
        return detect_sample_patterns(raw_df, sample_columns)
    
    def validate_sample_mapping(self, raw_df: pd.DataFrame, sample_condition_map: Dict[str, str]) -> Dict[str, List[str]]:
        """Validate sample condition mapping against raw data."""
        return validate_sample_mapping(raw_df, sample_condition_map)
    
    def summary(self) -> str:
        """Generate a summary of the payload configuration."""
        summary_lines = [
            "=== Payload Summary ===",
            f"Curtain Type: {self.curtain_type}",
            f"Search Groups: {len(self.get_search_groups())}",
            f"Conditions: {len(self.get_conditions())}",
            f"Samples: {len(self.get_samples())}"
        ]
        
        if self.get_search_groups():
            summary_lines.append("Search Groups:")
            for group in self.get_search_groups():
                color = self.data.get('settings', {}).get('colorMap', {}).get(group, 'N/A')
                summary_lines.append(f"  - {group} (Color: {color})")
        
        if self.get_conditions():
            summary_lines.append("Conditions:")
            for condition in self.get_conditions():
                color = self.data.get('settings', {}).get('colorMap', {}).get(condition, 'N/A')
                summary_lines.append(f"  - {condition} (Color: {color})")
        
        return "\n".join(summary_lines)
    
    # Factory methods
    @classmethod
    def from_dict(cls, data: Dict) -> 'Payload':
        """Create Payload instance from dictionary data."""
        curtain_type = data.get('curtain_type', 'TP')
        return cls(curtain_type=curtain_type, initial_payload=data)
    
    @classmethod  
    def from_json(cls, json_path: str) -> 'Payload':
        """Create Payload instance from JSON file."""
        data = parse_curtain_from_json(json_path)
        return cls.from_dict(data)
    
    @classmethod
    def create_with_data(cls, raw_df: pd.DataFrame, de_df: pd.DataFrame = None, 
                        curtain_type: str = 'TP', **kwargs) -> 'Payload':
        """
        Create Payload with automatic data processing and UniProt integration.
        
        Args:
            raw_df: Raw protein data DataFrame
            de_df: Differential expression DataFrame (required for PTM)
            curtain_type: Platform type ('TP' or 'PTM')
            **kwargs: Additional arguments for data processing
        """
        if curtain_type.upper() == 'PTM':
            if de_df is None:
                raise ValueError("de_df is required for CurtainPTM (curtain_type='PTM')")
            payload_dict = add_uniprot_data_ptm(deepcopy(curtain_base_payload), raw_df, de_df)
        else:
            payload_dict = add_uniprot_data(deepcopy(curtain_base_payload), raw_df)
        
        payload_dict['curtain_type'] = curtain_type.upper()
        return cls.from_dict(payload_dict)


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
    return json_data


def create_imputation_map(
    raw_df: pd.DataFrame | str, primary_id: str, sample_cols: List[str]
) -> Dict:
    """Create an imputation map indicating missing values in the dataset.
    
    This function identifies missing (null) values in sample columns for each
    primary identifier, creating a map that can be used to track which values
    were imputed during data processing.
    
    Args:
        raw_df: Raw data DataFrame or path to file containing raw data
        primary_id: Column name containing primary identifiers
        sample_cols: List of sample column names to check for missing values
        
    Returns:
        Dictionary mapping primary IDs to dictionaries of imputed columns.
        Structure: {primary_id: {column_name: True}} where True indicates
        the value was missing/imputed.
        
    Raises:
        ValueError: If specified columns are not found in the data
    """
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
    
    parser_columns = "accession,id,gene_names,protein_name,organism_name,organism_id,length,cc_subcellular_location,sequence,ft_var_seq,cc_alternative_products,ft_domain,xref_string,ft_mod_res,cc_function,cc_disease,cc_pharmaceutical,ft_mutagen,xref_mim"
    
    primary_id_raw = payload["rawForm"]["_primaryIDs"]
    primary_id_de = payload["differentialForm"]["_primaryIDs"]
    accession_col = payload["differentialForm"]["_accession"]
    if primary_id_raw not in raw_df.columns:
        raise ValueError(f"Primary ID column '{primary_id_raw}' not found in raw data")
    if primary_id_de not in de_df.columns:
        raise ValueError(f"Primary ID column '{primary_id_de}' not found in differential data")
    if accession_col not in de_df.columns:
        raise ValueError(f"Accession column '{accession_col}' not found in differential data")
    
    results = {}
    dataMap = {}
    db = {}
    organism = ""
    accMap = {}
    geneNameToPrimary = {}
    genesMap = {}
    
    accessionMap = {}
    accessionToPrimaryIDs = {}
    accessionList = []
    primaryIDsList = []
    
    acc_list = []
    uniprot_dataMap = {}
    
    primary_id_col_raw = raw_df[primary_id_raw].unique()
    primaryIDsList = list(primary_id_col_raw)
    
    accession_col_data = de_df[accession_col].dropna().unique()
    accessionList = list(accession_col_data)
    for accession in accessionList:
        if accession not in accessionMap:
            accessionMap[accession] = {}
            accessionMap[accession][accession] = True
        
        for split_acc in accession.split(";"):
            if split_acc not in accessionMap:
                accessionMap[split_acc] = {}
            accessionMap[split_acc][accession] = True
    
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
                
                dataMap[accession] = accession
                dataMap[primary_id] = accession
                if us.accession not in accMap:
                    accMap[us.accession] = [us.accession]
                else:
                    if us.accession not in accMap[us.accession]:
                        accMap[us.accession].append(us.accession)
                
                if us.accession not in uniprot_dataMap:
                    acc_list.append(us.accession)
    
    if acc_list:
        parser = UniprotParser(columns=parser_columns)
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
                        
                        if gene_name not in geneNameToPrimary:
                            geneNameToPrimary[gene_name] = {}
                        
                        if r["Entry"] in accessionToPrimaryIDs:
                            for primary_id in accessionToPrimaryIDs[r["Entry"]]:
                                geneNameToPrimary[gene_name][primary_id] = True
                    
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
    
    payload["extraData"] = {
        "uniprot": {
            "accMap": {"dataType": "Map", "value": list(accMap.items())},
            "dataMap": {"dataType": "Map", "value": list(uniprot_dataMap.items())},
            "db": {"dataType": "Map", "value": list(db.items())},
            "organism": organism,
            "geneNameToPrimary": geneNameToPrimary,
            "results": {"dataType": "Map", "value": list(results.items())},
        },
        "data": {
            "accessionToPrimaryIDs": accessionToPrimaryIDs,
            "primaryIDsList": primaryIDsList,
            "accessionList": accessionList,
            "accessionMap": accessionMap,
            "genesMap": genesMap,
            "allGenes": allGenes,
            "dataMap": {"dataType": "Map", "value": list(dataMap.items())},
        },
    }
    
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
    
    if 'volcanoAxis' not in payload['settings']:
        payload['settings']['volcanoAxis'] = {
            'minX': None, 'maxX': None, 'minY': None, 'maxY': None,
            'x': "Log2FC", 'y': "-log10(p-value)",
            'dtickX': None, 'dtickY': None,
            'ticklenX': 5, 'ticklenY': 5
        }
    
    if 'x_min' in kwargs:
        payload['settings']['volcanoAxis']['minX'] = kwargs['x_min']
    if 'x_max' in kwargs:
        payload['settings']['volcanoAxis']['maxX'] = kwargs['x_max']
    if 'y_min' in kwargs:
        payload['settings']['volcanoAxis']['minY'] = kwargs['y_min']
    if 'y_max' in kwargs:
        payload['settings']['volcanoAxis']['maxY'] = kwargs['y_max']
    
    if 'x_title' in kwargs:
        payload['settings']['volcanoAxis']['x'] = kwargs['x_title']
    if 'y_title' in kwargs:
        payload['settings']['volcanoAxis']['y'] = kwargs['y_title']
    
    if 'x_tick_interval' in kwargs:
        payload['settings']['volcanoAxis']['dtickX'] = kwargs['x_tick_interval']
    if 'y_tick_interval' in kwargs:
        payload['settings']['volcanoAxis']['dtickY'] = kwargs['y_tick_interval']
    if 'x_tick_length' in kwargs:
        payload['settings']['volcanoAxis']['ticklenX'] = kwargs['x_tick_length']
    if 'y_tick_length' in kwargs:
        payload['settings']['volcanoAxis']['ticklenY'] = kwargs['y_tick_length']
    
    if 'volcanoPlotDimension' not in payload['settings']:
        payload['settings']['volcanoPlotDimension'] = {
            'width': 800, 'height': 1000,
            'margin': {'l': None, 'r': None, 'b': None, 't': None}
        }
    
    if 'width' in kwargs:
        payload['settings']['volcanoPlotDimension']['width'] = kwargs['width']
    if 'height' in kwargs:
        payload['settings']['volcanoPlotDimension']['height'] = kwargs['height']
    
    if any(k in kwargs for k in ['margin_left', 'margin_right', 'margin_bottom', 'margin_top']):
        if 'margin_left' in kwargs:
            payload['settings']['volcanoPlotDimension']['margin']['l'] = kwargs['margin_left']
        if 'margin_right' in kwargs:
            payload['settings']['volcanoPlotDimension']['margin']['r'] = kwargs['margin_right']
        if 'margin_bottom' in kwargs:
            payload['settings']['volcanoPlotDimension']['margin']['b'] = kwargs['margin_bottom']
        if 'margin_top' in kwargs:
            payload['settings']['volcanoPlotDimension']['margin']['t'] = kwargs['margin_top']
    
    if 'volcanoPlotGrid' not in payload['settings']:
        payload['settings']['volcanoPlotGrid'] = {'x': True, 'y': True}
    
    if 'show_x_grid' in kwargs:
        payload['settings']['volcanoPlotGrid']['x'] = kwargs['show_x_grid']
    if 'show_y_grid' in kwargs:
        payload['settings']['volcanoPlotGrid']['y'] = kwargs['show_y_grid']
    
    if 'title' in kwargs:
        payload['settings']['volcanoPlotTitle'] = kwargs['title']
    if 'legend_x' in kwargs:
        payload['settings']['volcanoPlotLegendX'] = kwargs['legend_x']
    if 'legend_y' in kwargs:
        payload['settings']['volcanoPlotLegendY'] = kwargs['legend_y']
    
    if 'marker_size' in kwargs:
        payload['settings']['scatterPlotMarkerSize'] = kwargs['marker_size']
    
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
    
    if 'columnSize' not in payload['settings']:
        payload['settings']['columnSize'] = {
            'barChart': 0, 'averageBarChart': 0,
            'violinPlot': 0, 'profilePlot': 0
        }
    
    if 'bar_chart_width' in kwargs:
        payload['settings']['columnSize']['barChart'] = kwargs['bar_chart_width']
    if 'average_bar_chart_width' in kwargs:
        payload['settings']['columnSize']['averageBarChart'] = kwargs['average_bar_chart_width']
    if 'violin_plot_width' in kwargs:
        payload['settings']['columnSize']['violinPlot'] = kwargs['violin_plot_width']
    if 'profile_plot_width' in kwargs:
        payload['settings']['columnSize']['profilePlot'] = kwargs['profile_plot_width']
    
    if 'condition_colors' in kwargs:
        if 'barchartColorMap' not in payload['settings']:
            payload['settings']['barchartColorMap'] = {}
        payload['settings']['barchartColorMap'].update(kwargs['condition_colors'])
    
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
    
    if 'font_family' in kwargs:
        payload['settings']['plotFontFamily'] = kwargs['font_family']
    
    if 'p_cutoff' in kwargs:
        payload['settings']['pCutoff'] = kwargs['p_cutoff']
    if 'fc_cutoff' in kwargs:
        payload['settings']['log2FCCutoff'] = kwargs['fc_cutoff']
    
    if 'default_colors' in kwargs:
        payload['settings']['defaultColorList'] = kwargs['default_colors']
    
    if 'condition_colors' in kwargs:
        if 'colorMap' not in payload['settings']:
            payload['settings']['colorMap'] = {}
        payload['settings']['colorMap'].update(kwargs['condition_colors'])
    
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
    
    if 'custom_ptm_data' in kwargs:
        payload['settings']['customPTMData'] = kwargs['custom_ptm_data']
    
    if 'variant_corrections' in kwargs:
        payload['settings']['variantCorrection'] = kwargs['variant_corrections']
    
    if 'custom_sequences' in kwargs:
        payload['settings']['customSequences'] = kwargs['custom_sequences']
    
    return payload


def configure_sample_conditions(payload: Dict, sample_condition_map: Dict[str, str], 
                              condition_colors: Dict[str, str] = None, 
                              condition_order: List[str] = None) -> Dict:
    """
    Configure custom sample-to-condition mapping, overriding the default condition.replicate parsing.
    
    This function allows you to assign samples to custom condition groups instead of relying on
    the default naming convention (condition.replicate format).
    
    Args:
        payload (Dict): The payload to modify.
        sample_condition_map (Dict[str, str]): Dictionary mapping sample names to condition names.
                                             Keys are sample names, values are condition names.
        condition_colors (Dict[str, str], optional): Dictionary mapping condition names to hex colors.
        condition_order (List[str], optional): List specifying the order of conditions in visualizations.
    
    Returns:
        Dict: Updated payload with custom sample-condition mapping.
    
    Example:
        >>> # Custom sample assignment
        >>> sample_mapping = {
        >>>     'Sample_001': 'Control',
        >>>     'Sample_002': 'Control', 
        >>>     'Sample_003': 'Treatment_A',
        >>>     'Sample_004': 'Treatment_A',
        >>>     'Sample_005': 'Treatment_B',
        >>>     'Sample_006': 'Treatment_B'
        >>> }
        >>> 
        >>> condition_colors = {
        >>>     'Control': '#808080',
        >>>     'Treatment_A': '#FF6B6B',
        >>>     'Treatment_B': '#4ECDC4'
        >>> }
        >>>
        >>> payload = configure_sample_conditions(payload, sample_mapping, condition_colors)
    """
    if 'settings' not in payload:
        payload['settings'] = {}
    
    # Initialize sample mapping structures
    sample_map = {}
    sample_order = {}
    conditions = []
    
    # Create sample map entries and collect conditions
    for sample_name, condition_name in sample_condition_map.items():
        if condition_name not in conditions:
            conditions.append(condition_name)
        
        # Extract or generate replicate number
        if '.' in sample_name:
            # If sample follows condition.replicate format, extract replicate
            parts = sample_name.split('.')
            replicate = parts[-1]
        else:
            # Generate sequential replicate numbers for samples in each condition
            existing_samples = [s for s, c in sample_condition_map.items() if c == condition_name]
            replicate = str(existing_samples.index(sample_name) + 1)
        
        # Create sample map entry
        sample_map[sample_name] = {
            'replicate': replicate,
            'condition': condition_name,
            'name': sample_name
        }
        
        # Add to sample order
        if condition_name not in sample_order:
            sample_order[condition_name] = []
        if sample_name not in sample_order[condition_name]:
            sample_order[condition_name].append(sample_name)
    
    # Update payload settings
    payload['settings']['sampleMap'] = sample_map
    payload['settings']['sampleOrder'] = sample_order
    
    # Set condition order
    if condition_order is not None:
        payload['settings']['conditionOrder'] = condition_order
    else:
        payload['settings']['conditionOrder'] = conditions
    
    # Set sample visibility (all visible by default)
    if 'sampleVisible' not in payload['settings']:
        payload['settings']['sampleVisible'] = {}
    for sample_name in sample_condition_map.keys():
        if sample_name not in payload['settings']['sampleVisible']:
            payload['settings']['sampleVisible'][sample_name] = True
    
    # Set condition colors
    if 'colorMap' not in payload['settings']:
        payload['settings']['colorMap'] = {}
    
    if condition_colors:
        payload['settings']['colorMap'].update(condition_colors)
    else:
        # Auto-assign colors using default color list
        default_colors = [
            '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', 
            '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf'
        ]
        for i, condition in enumerate(conditions):
            if condition not in payload['settings']['colorMap']:
                payload['settings']['colorMap'][condition] = default_colors[i % len(default_colors)]
    
    # Set bar chart colors (copy from main color map)
    if 'barchartColorMap' not in payload['settings']:
        payload['settings']['barchartColorMap'] = {}
    payload['settings']['barchartColorMap'].update(payload['settings']['colorMap'])
    
    return payload


def detect_sample_patterns(raw_df: pd.DataFrame, sample_columns: List[str]) -> Dict[str, str]:
    """
    Automatically detect potential condition patterns from sample column names.
    
    This function analyzes sample names and suggests condition mappings based on common patterns:
    - Standard condition.replicate format (e.g., "Control.1", "Treatment.2")  
    - Underscore patterns (e.g., "Control_rep1", "Treatment_A_rep2")
    - Prefix patterns (e.g., "Ctrl01", "Trt01", "Drug01")
    - Numeric patterns (e.g., "Sample1", "Sample2" -> groups by ranges)
    
    Args:
        raw_df (pd.DataFrame): Raw data DataFrame containing sample columns.
        sample_columns (List[str]): List of sample column names to analyze.
    
    Returns:
        Dict[str, str]: Suggested mapping of sample names to condition names.
        
    Example:
        >>> raw_df = pd.read_csv("data.csv")
        >>> samples = ["Control.1", "Control.2", "Treatment.1", "Treatment.2"] 
        >>> suggestions = detect_sample_patterns(raw_df, samples)
        >>> print(suggestions)
        {'Control.1': 'Control', 'Control.2': 'Control', 
         'Treatment.1': 'Treatment', 'Treatment.2': 'Treatment'}
    """
    suggestions = {}
    
    for sample in sample_columns:
        condition = None
        
        # Pattern 1: Standard condition.replicate format
        if '.' in sample:
            parts = sample.split('.')
            if len(parts) >= 2 and parts[-1].isdigit():
                condition = '.'.join(parts[:-1])
        
        # Pattern 2: Underscore patterns (Control_rep1, Treatment_A_rep2)
        elif '_rep' in sample.lower() or '_r' in sample.lower():
            if '_rep' in sample.lower():
                condition = sample.lower().split('_rep')[0]
            elif '_r' in sample.lower():
                condition = sample.lower().split('_r')[0]
            condition = condition.replace('_', ' ').title()
        
        # Pattern 3: Common prefix patterns
        elif any(sample.lower().startswith(prefix) for prefix in ['ctrl', 'control', 'con']):
            condition = 'Control'
        elif any(sample.lower().startswith(prefix) for prefix in ['trt', 'treat', 'treatment']):
            condition = 'Treatment'  
        elif any(sample.lower().startswith(prefix) for prefix in ['drug', 'compound']):
            condition = 'Drug'
        elif sample.lower().startswith('mock'):
            condition = 'Mock'
        elif sample.lower().startswith('vehicle'):
            condition = 'Vehicle'
        
        # Pattern 4: Numeric grouping (Sample1-3 = Group1, Sample4-6 = Group2, etc.)
        elif sample.lower().startswith('sample') and sample[-1].isdigit():
            sample_num = int(sample.split('sample')[1] if 'sample' in sample.lower() else sample[-1])
            group_num = (sample_num - 1) // 3 + 1  # Group every 3 samples
            condition = f'Group_{group_num}'
        
        # Pattern 5: Generic alphanumeric - try to extract meaningful parts
        else:
            # Remove common suffixes/numbers to find base condition
            base = sample
            for suffix in ['1', '2', '3', '01', '02', '03', '_1', '_2', '_3']:
                if base.endswith(suffix):
                    base = base[:-len(suffix)]
                    break
            condition = base if base else f'Condition_{sample}'
        
        suggestions[sample] = condition if condition else 'Unknown'
    
    return suggestions


def validate_sample_mapping(raw_df: pd.DataFrame, sample_condition_map: Dict[str, str]) -> Dict[str, List[str]]:
    """
    Validate a sample-condition mapping against the raw data.
    
    Checks for common issues:
    - Missing samples (in mapping but not in data)
    - Extra samples (in data but not in mapping) 
    - Unbalanced conditions (significantly different sample counts)
    - Invalid sample names
    
    Args:
        raw_df (pd.DataFrame): Raw data DataFrame.
        sample_condition_map (Dict[str, str]): Sample to condition mapping to validate.
    
    Returns:
        Dict[str, List[str]]: Dictionary with validation results:
            - 'missing_samples': Samples in mapping but not in data
            - 'extra_samples': Samples in data but not in mapping  
            - 'invalid_samples': Samples with invalid names/characters
            - 'unbalanced_conditions': Conditions with very different sample counts
            - 'warnings': General warnings about the mapping
    """
    validation_results = {
        'missing_samples': [],
        'extra_samples': [], 
        'invalid_samples': [],
        'unbalanced_conditions': [],
        'warnings': []
    }
    
    data_columns = raw_df.columns.tolist()
    mapping_samples = set(sample_condition_map.keys())
    data_samples = set(data_columns)
    
    # Check for missing and extra samples
    validation_results['missing_samples'] = list(mapping_samples - data_samples)
    validation_results['extra_samples'] = list(data_samples - mapping_samples)
    
    # Check for invalid sample names
    for sample in mapping_samples:
        if not sample.replace('_', '').replace('.', '').replace('-', '').isalnum():
            validation_results['invalid_samples'].append(sample)
    
    # Check condition balance
    condition_counts = {}
    for condition in sample_condition_map.values():
        condition_counts[condition] = condition_counts.get(condition, 0) + 1
    
    if len(condition_counts) > 1:
        counts = list(condition_counts.values())
        max_count, min_count = max(counts), min(counts)
        if max_count > min_count * 2:  # Flag if 2x+ difference
            validation_results['unbalanced_conditions'] = [
                f"{condition}: {count} samples" for condition, count in condition_counts.items()
            ]
    
    # Generate warnings
    if validation_results['missing_samples']:
        validation_results['warnings'].append(f"Missing {len(validation_results['missing_samples'])} samples from data")
    
    if validation_results['extra_samples']:
        validation_results['warnings'].append(f"Found {len(validation_results['extra_samples'])} unmapped samples in data")
    
    if validation_results['unbalanced_conditions']:
        validation_results['warnings'].append("Conditions have unbalanced sample counts")
    
    return validation_results


def create_batch_search_group(payload: Dict, group_name: str, search_input: str, 
                             search_type: str = "Gene Names", color: str = None,
                             specific_comparison_label: str = None) -> Dict:
    """
    Replicate exact frontend batch search logic from batch-search.component + protein-selections.component.
    
    Args:
        search_input: Multi-line string with protein/gene identifiers (like frontend textarea)
        search_type: "Gene Names" or "Primary IDs" (TP) / "Accession IDs" (PTM)
    """
    print(f"DEBUG BATCH: Starting batch search for '{group_name}'")
    print(f"DEBUG BATCH: Search type: {search_type}")
    print(f"DEBUG BATCH: Input: {repr(search_input)}")
    
    is_curtain_ptm = payload.get('curtain_type') == 'PTM'
    print(f"DEBUG BATCH: Is CurtainPTM: {is_curtain_ptm}")
    
    # Step 1: Process input like frontend handleSubmit()
    batch_data = {}
    for line in search_input.replace("\r", "").split("\n"):
        line = line.strip().upper()
        if line:
            parts = line.split(";")
            batch_data[line] = []
            for part in parts:
                part = part.strip()
                if part:
                    batch_data[line].append(part)
    
    print(f"DEBUG BATCH: Processed batch_data: {batch_data}")
    
    # Step 2: Process like frontend getPrimaryIDsDataFromBatch()
    result_primary_ids = []
    
    for input_line in batch_data:
        print(f"DEBUG BATCH: Processing line: {input_line}")
        
        # Try exact match first
        primary_ids = _parse_data_frontend_style(payload, input_line, search_type, exact=True)
        print(f"DEBUG BATCH: Exact match result: {primary_ids}")
        
        # If no exact match, try inexact matching with individual parts
        if not primary_ids:
            for input_part in batch_data[input_line]:
                print(f"DEBUG BATCH: Trying inexact match for part: {input_part}")
                primary_ids = _parse_data_frontend_style(payload, input_part, search_type, exact=False)
                print(f"DEBUG BATCH: Inexact match result: {primary_ids}")
                if primary_ids:
                    break
        
        # Add found primary IDs
        result_primary_ids.extend(primary_ids)
    
    print(f"DEBUG BATCH: Final result_primary_ids: {result_primary_ids}")
    
    # Step 3: Create search group with found primary IDs
    if result_primary_ids:
        return add_search_group(payload, group_name, protein_ids=result_primary_ids, 
                              color=color, specific_comparison_label=specific_comparison_label)
    else:
        print(f"No proteins found for batch search: {group_name}")
        return payload

def _parse_data_frontend_style(payload: Dict, identifier: str, search_type: str, exact: bool) -> List[str]:
    """Replicate frontend parseData() logic exactly."""
    print(f"DEBUG PARSE: identifier='{identifier}', search_type='{search_type}', exact={exact}")
    
    is_curtain_ptm = payload.get('curtain_type') == 'PTM'
    
    if search_type == "Gene Names":
        if exact:
            print(f"DEBUG PARSE: Doing exact gene name lookup")
            # Direct gene name lookup
            if is_curtain_ptm:
                result = _get_primary_from_gene_names_ptm(payload, identifier)
                print(f"DEBUG PARSE: PTM gene lookup result: {result}")
                return result
            else:
                result = _get_primary_ids_from_gene_names_tp(payload, identifier)
                print(f"DEBUG PARSE: TP gene lookup result: {result}")
                return result
        else:
            print(f"DEBUG PARSE: Doing fuzzy gene name lookup")
            # Fuzzy lookup using genesMap
            genes_map = payload.get('extraData', {}).get('data', {}).get('genesMap', {})
            print(f"DEBUG PARSE: genesMap has {len(genes_map)} entries")
            if identifier in genes_map:
                gene_dict = genes_map[identifier]
                if isinstance(gene_dict, dict):
                    for mapped_gene in gene_dict.keys():
                        if is_curtain_ptm:
                            result = _get_primary_from_gene_names_ptm(payload, mapped_gene)
                        else:
                            result = _get_primary_ids_from_gene_names_tp(payload, mapped_gene)
                        if result:
                            return result
    
    elif search_type in ["Primary IDs", "Accession IDs"]:
        print(f"DEBUG PARSE: Processing protein/accession ID")
        if exact:
            # Direct protein ID lookup
            if is_curtain_ptm:
                return _get_primary_from_acc_ptm(payload, identifier)
            else:
                return _get_primary_ids_from_acc_tp(payload, identifier)
        else:
            # Fuzzy lookup using primaryIDsMap/accessionMap
            if is_curtain_ptm:
                acc_map = payload.get('extraData', {}).get('data', {}).get('accessionMap', {})
                if identifier in acc_map:
                    acc_dict = acc_map[identifier]
                    if isinstance(acc_dict, dict):
                        for mapped_acc in acc_dict.keys():
                            result = _get_primary_from_acc_ptm(payload, mapped_acc)
                            if result:
                                return result
            else:
                primary_ids_map = payload.get('extraData', {}).get('data', {}).get('primaryIDsmap', {})
                if identifier in primary_ids_map:
                    primary_dict = primary_ids_map[identifier]
                    if isinstance(primary_dict, dict):
                        for mapped_primary in primary_dict.keys():
                            result = _get_primary_ids_from_acc_tp(payload, mapped_primary)
                            if result:
                                return result
    
    print(f"DEBUG PARSE: No results found, returning empty list")
    return []

def add_search_group_with_type(payload: Dict, group_name: str, curtain_type: str, protein_ids: List[str] = None,
                    gene_names: List[str] = None, color: str = None, 
                    specific_comparison_label: str = None, raw_df: pd.DataFrame = None,
                    primary_id_col: str = None, gene_name_col: str = None) -> Dict:
    """
    Add a search group to the payload following frontend comparison label logic.
    
    Behavior matches frontend exactly:
    - Curtain: Creates groups for ALL selected comparison labels with " (label)" suffix
    - CurtainPTM: Creates single group without comparison suffix
    - Default: Uses "CurtainSetComparison" column with label "1" if no comparison set
    
    Args:
        payload: The payload to modify.
        group_name: Name of the search group.
        curtain_type: Platform type ('TP' for Curtain, 'PTM' for CurtainPTM).
        protein_ids: List of protein/accession IDs to include in the group.
        gene_names: List of gene names to include in the group.
        color: Hex color code for the group. Auto-assigns if not provided.
        specific_comparison_label: Override to create group for specific label only.
        raw_df: Raw data for fallback searching when using gene names.
        primary_id_col: Primary ID column name for raw data searching.
        gene_name_col: Gene name column for raw data searching.
    
    Returns:
        Updated payload with the new search group(s).
    """
    if 'selected' not in payload:
        payload['selected'] = []
    if 'selectedMap' not in payload:
        payload['selectedMap'] = {}
    if 'selectOperationNames' not in payload:
        payload['selectOperationNames'] = []
    if 'settings' not in payload:
        payload['settings'] = {}
    if 'colorMap' not in payload['settings']:
        payload['settings']['colorMap'] = {}
    
    is_curtain_ptm = curtain_type.upper() == 'PTM'
    
    # Get comparison configuration
    diff_form = payload.get('differentialForm', {})
    comparison_col = diff_form.get('_comparison', '')
    comparison_select = diff_form.get('_comparisonSelect', [])
    
    # Determine which labels to use for group creation
    labels_to_use = []
    
    if specific_comparison_label:
        # User specified a specific label
        labels_to_use = [specific_comparison_label]
    elif is_curtain_ptm:
        # CurtainPTM: No comparison suffix, single group
        labels_to_use = []
    elif comparison_col == 'CurtainSetComparison' or not comparison_col:
        # Default case: use "1" as label
        labels_to_use = ['1']
    elif isinstance(comparison_select, list) and comparison_select:
        # Curtain with multiple selected comparisons
        labels_to_use = comparison_select
    elif isinstance(comparison_select, str) and comparison_select:
        # Single comparison selected
        labels_to_use = [comparison_select]
    else:
        # Fallback: use "1" as default
        labels_to_use = ['1']
    
    # Validate input - at least one of protein_ids or gene_names must be provided
    if not protein_ids and not gene_names:
        raise ValueError("Either protein_ids or gene_names must be provided")
    
    # Build the final list of protein IDs
    all_protein_ids = []
    
    # Add direct protein IDs
    if protein_ids:
        matched_proteins = _search_proteins_frontend_style(payload, protein_ids, raw_df, primary_id_col)
        all_protein_ids.extend(matched_proteins)
    
    if gene_names:
        matched_proteins = _search_genes_frontend_style(payload, gene_names, is_curtain_ptm, raw_df, primary_id_col, gene_name_col)
        all_protein_ids.extend(matched_proteins)
    
    final_protein_ids = list(dict.fromkeys(all_protein_ids))
    
    if not final_protein_ids:
        return payload
    
    if is_curtain_ptm:
        final_group_name = group_name
        _add_single_search_group(payload, final_group_name, final_protein_ids, color)
    else:
        for label in labels_to_use:
            final_group_name = f"{group_name} ({label})"
            _add_single_search_group(payload, final_group_name, final_protein_ids, color)
    
    return payload


def add_search_group(payload: Dict, group_name: str, curtain_type: str, protein_ids: List[str] = None,
                    gene_names: List[str] = None, color: str = None, 
                    specific_comparison_label: str = None, raw_df: pd.DataFrame = None,
                    primary_id_col: str = None, gene_name_col: str = None) -> Dict:
    """
    Add a search group to the payload.
    
    Args:
        payload: The payload to modify.
        group_name: Name of the search group.
        curtain_type: Platform type ('TP' for Curtain, 'PTM' for CurtainPTM).
        protein_ids: List of protein/accession IDs to include in the group.
        gene_names: List of gene names to include in the group.
        color: Hex color code for the group. Auto-assigns if not provided.
        specific_comparison_label: Override to create group for specific label only.
        raw_df: Raw data for fallback searching when using gene names.
        primary_id_col: Primary ID column name for raw data searching.
        gene_name_col: Gene name column for raw data searching.
    
    Returns:
        Updated payload with the new search group(s).
    """
    return add_search_group_with_type(payload, group_name, curtain_type, protein_ids, gene_names, color,
                                    specific_comparison_label, raw_df, primary_id_col, gene_name_col)


def _add_single_search_group(payload: Dict, group_name: str, primary_ids: List[str], color: str = None):
    # Add group name to operation list if not already present
    if group_name not in payload['selectionsName']:
        payload['selectionsName'].append(group_name)
    
    # Assign color from default palette if not provided
    if color is None:
        default_colors = [
            "#fd7f6f", "#7eb0d5", "#b2e061", "#bd7ebe", 
            "#ffb55a", "#ffee65", "#beb9db", "#fdcce5", "#8bd3c7"
        ]
        
        # Find next available color
        used_colors = set(payload['settings']['colorMap'].values())
        available_colors = [c for c in default_colors if c not in used_colors]
        
        if available_colors:
            color = available_colors[0]
        else:
            # If all default colors used, cycle through them
            color_index = len(payload['settings']['colorMap']) % len(default_colors)
            color = default_colors[color_index]
    
    # Assign color to the group
    payload['settings']['colorMap'][group_name] = color
    
    for primary_id in primary_ids:
        if primary_id not in payload['selections']:
            payload['selections'].append(primary_id)
    
    for primary_id in primary_ids:
        if primary_id not in payload['selectionsMap']:
            payload['selectionsMap'][primary_id] = {}
        payload['selectionsMap'][primary_id][group_name] = True


def add_multiple_search_groups(payload: Dict, groups: Dict[str, Dict]) -> Dict:
    """
    Add multiple search groups to the payload at once.
    
    Args:
        payload (Dict): The payload to modify.
        groups (Dict[str, Dict]): Dictionary where keys are group names and values contain:
            - 'proteins' (List[str]): List of protein IDs
            - 'color' (str, optional): Hex color code
            - 'specific_comparison_label' (str, optional): Specific comparison label override
    
    Returns:
        Dict: Updated payload with all search groups.
    
    Example:
        >>> groups = {
        ...     "Kinases": {
        ...         "proteins": ["P12345", "Q67890"],
        ...         "color": "#FF6B6B"
        ...     },
        ...     "Phosphatases": {
        ...         "proteins": ["P11111", "Q22222"],
        ...         "color": "#4ECDC4",
        ...         "comparison": "Treatment_vs_Control"
        ...     }
        ... }
        >>> payload = add_multiple_search_groups(payload, groups)
    """
    for group_name, group_config in groups.items():
        # Support both 'proteins'/'protein_ids' and 'gene_names' keys
        protein_ids = group_config.get('proteins') or group_config.get('protein_ids')
        gene_names = group_config.get('gene_names')
        
        payload = add_search_group(
            payload,
            group_name,
            protein_ids=protein_ids,
            gene_names=gene_names,
            color=group_config.get('color'),
            specific_comparison_label=group_config.get('specific_comparison_label')
        )
    
    return payload


def create_search_group_from_gene_list(payload: Dict, group_name: str, gene_names: List[str],
                                     raw_df: pd.DataFrame = None, primary_id_col: str = None,
                                     gene_name_col: str = None, color: str = None,
                                     specific_comparison_label: str = None) -> Dict:
    """
    Create a search group by matching gene names using the exact frontend search logic.
    
    This function replicates the frontend batch search functionality:
    - Converts all gene names to uppercase for matching
    - Handles semicolon-separated gene names  
    - Uses exact matching through genesMap and geneNameToAcc/geneNameToPrimary
    - Falls back to raw data searching if no UniProt data available
    
    Args:
        payload (Dict): The payload to modify.
        group_name (str): Name of the search group.
        gene_names (List[str]): List of gene names to search for (case-insensitive).
        raw_df (pd.DataFrame, optional): Raw data DataFrame. Uses payload data if not provided.
        primary_id_col (str, optional): Primary ID column. Uses payload config if not provided.
        gene_name_col (str, optional): Gene name column. Uses payload config if not provided.
        color (str, optional): Hex color code for the group.
        specific_comparison_label (str, optional): Specific comparison label override.
    
    Returns:
        Dict: Updated payload with the new search group containing matched proteins.
        
    Example:
        >>> # Using frontend-style gene search
        >>> payload = create_search_group_from_gene_list(
        ...     payload, 
        ...     "Cancer Genes",
        ...     ["BRCA1", "TP53", "MYC;EGFR"]  # Supports semicolon separation
        ... )
    """
    matched_proteins = _search_genes_frontend_style(payload, gene_names, raw_df, primary_id_col, gene_name_col)
    
    if not matched_proteins:
        print(f"Warning: No proteins found for genes {gene_names}")
        return payload
    
    return add_search_group(payload, group_name, matched_proteins, color, specific_comparison_label)


def create_search_group_from_protein_list(payload: Dict, group_name: str, protein_ids: List[str],
                                        raw_df: pd.DataFrame = None, primary_id_col: str = None,
                                        color: str = None, specific_comparison_label: str = None) -> Dict:
    """
    Create a search group by matching protein IDs using the exact frontend search logic.
    
    This function replicates the frontend protein ID search functionality:
    - Handles semicolon-separated protein IDs
    - Uses exact matching through primaryIDsMap and accessionToPrimaryIDs
    - Supports both primary IDs and accession IDs
    
    Args:
        payload (Dict): The payload to modify.
        group_name (str): Name of the search group.
        protein_ids (List[str]): List of protein IDs to search for.
        raw_df (pd.DataFrame, optional): Raw data DataFrame. Uses payload data if not provided.
        primary_id_col (str, optional): Primary ID column. Uses payload config if not provided.
        color (str, optional): Hex color code for the group.
        specific_comparison_label (str, optional): Specific comparison label override.
    
    Returns:
        Dict: Updated payload with the new search group containing matched proteins.
        
    Example:
        >>> # Using frontend-style protein ID search  
        >>> payload = create_search_group_from_protein_list(
        ...     payload,
        ...     "Target Proteins", 
        ...     ["P12345", "Q67890;P11111"]  # Supports semicolon separation
        ... )
    """
    matched_proteins = _search_proteins_frontend_style(payload, protein_ids, raw_df, primary_id_col)
    
    if not matched_proteins:
        print(f"Warning: No proteins found for IDs {protein_ids}")
        return payload
    
    return add_search_group(payload, group_name, matched_proteins, color, specific_comparison_label)


def _search_genes_frontend_style(payload: Dict, gene_names: List[str], is_curtain_ptm: bool, raw_df: pd.DataFrame = None,
                                primary_id_col: str = None, gene_name_col: str = None) -> List[str]:
    """
    Convert gene names to primary IDs using frontend logic.
    
    Args:
        payload: The payload containing gene mapping data.
        gene_names: List of gene names to search for.
        is_curtain_ptm: Whether this is CurtainPTM (True) or regular Curtain (False).
        raw_df: Raw data for fallback searching.
        primary_id_col: Primary ID column name for raw data searching.
        gene_name_col: Gene name column for raw data searching.
    
    Returns:
        List of primary IDs found for the gene names.
    """
    result = []
    
    for gene_input in gene_names:
        for gene_part in gene_input.split(';'):
            gene_name = gene_part.strip().upper()
            if not gene_name:
                continue
                
            if is_curtain_ptm:
                primary_ids = _get_primary_from_gene_names_ptm(payload, gene_name)
            else:
                primary_ids = _get_primary_ids_from_gene_names_tp(payload, gene_name)
            
            for primary_id in primary_ids:
                if primary_id not in result:
                    result.append(primary_id)
    
    return result

def _get_primary_from_gene_names_ptm(payload: Dict, gene_name: str) -> List[str]:
    """
    CurtainPTM gene name to primary ID conversion.
    
    Args:
        payload: The payload containing gene mapping data.
        gene_name: Gene name to search for.
    
    Returns:
        List of primary IDs found for the gene name.
    """
    result = []
    if 'extraData' not in payload or 'uniprot' not in payload['extraData']:
        return result
    
    gene_name_to_primary = payload['extraData']['uniprot'].get('geneNameToPrimary', {})
    gene_name_upper = gene_name.upper()
    
    if gene_name_upper in gene_name_to_primary:
        primary_dict = gene_name_to_primary[gene_name_upper]
        
        if isinstance(primary_dict, dict):
            for primary_id in primary_dict.keys():
                if primary_id not in result:
                    result.append(primary_id)
    
    return result

def _get_primary_ids_from_gene_names_tp(payload: Dict, gene_name: str) -> List[str]:
    """
    Curtain TP gene name to primary ID conversion.
    Logic: gene_name -> geneNameToAcc -> accessions -> primaryIDsMap -> final results
    
    Args:
        payload: The payload containing gene mapping data.
        gene_name: Gene name to search for.
    
    Returns:
        List of primary IDs found for the gene name.
    """
    result = []
    
    if 'extraData' not in payload or 'uniprot' not in payload['extraData']:
        return result
    
    uniprot_data = payload['extraData']['uniprot']
    gene_name_to_acc = uniprot_data.get('geneNameToAcc', {})
    
    data_section = payload.get('extraData', {}).get('data', {})
    primary_ids_map = data_section.get('primaryIDsMap', {})  
    
    gene_name_upper = gene_name.upper()
    
    if gene_name_upper in gene_name_to_acc:
        acc_dict = gene_name_to_acc[gene_name_upper]
        
        if isinstance(acc_dict, dict):
            for acc in acc_dict.keys():
                if acc in primary_ids_map:
                    primary_dict = primary_ids_map[acc]
                    if isinstance(primary_dict, dict):
                        for primary_id in primary_dict.keys():
                            if primary_id not in result:
                                result.append(primary_id)
    
    return result


def create_significance_groups(payload: Dict, curtain_type: str, p_cutoff: float = 0.05, fc_cutoff: float = 0.6, 
                              colors: Dict[str, str] = None) -> Dict:
    """
    Create automatic significance groups based on p-value and fold change cutoffs.
    
    This replicates the frontend's automatic grouping system that categorizes
    datapoints into 4 groups based on significance and fold change thresholds.
    
    Args:
        payload: The payload to modify.
        curtain_type: Platform type ('TP' for Curtain, 'PTM' for CurtainPTM).
        p_cutoff: P-value threshold for significance (default: 0.05).
        fc_cutoff: Log2 fold change threshold (default: 0.6).
        colors: Optional dictionary mapping group patterns to colors.
               Keys should be like "P-value > FC > " format.
    
    Returns:
        Updated payload with significance group settings.
    """
    if 'settings' not in payload:
        payload['settings'] = {}
    
    payload['settings']['pCutoff'] = p_cutoff
    payload['settings']['log2FCCutoff'] = fc_cutoff
    
    if colors:
        if 'colorMap' not in payload['settings']:
            payload['settings']['colorMap'] = {}
        
        _apply_significance_group_colors(payload, colors, p_cutoff, fc_cutoff, curtain_type)
    
    return payload


def get_significance_group(x: float, y: float, p_cutoff: float = 0.05, fc_cutoff: float = 0.6, 
                          comparison_label: str = "1", curtain_type: str = "TP") -> tuple:
    """
    Determine significance group for a datapoint based on x,y coordinates.
    
    Replicates the frontend significantGroup() function logic with comparison labels.
    
    Args:
        x: Log2 fold change value.
        y: -log10(p-value) value.
        p_cutoff: P-value threshold for significance.
        fc_cutoff: Log2 fold change threshold.
        comparison_label: Comparison label to append to group name.
        curtain_type: Platform type ('TP' or 'PTM') for group naming.
    
    Returns:
        Tuple of (full_group_name, position_key) where:
        - full_group_name: Complete group name with comparison label suffix
        - position_key: Position identifier for consistent coloring
    """
    import math
    ylog = -math.log10(p_cutoff)
    groups = []
    position = ""
    
    if ylog > y:
        groups.append(f"P-value > {p_cutoff}")
        position = "P-value > "
    else:
        groups.append(f"P-value <= {p_cutoff}")
        position = "P-value <= "
    
    if abs(x) > fc_cutoff:
        groups.append(f"FC > {fc_cutoff}")
        position += "FC > "
    else:
        groups.append(f"FC <= {fc_cutoff}")
        position += "FC <= "
    
    base_group_name = ";".join(groups)
    
    # Add comparison suffix based on curtain type
    if curtain_type.upper() == 'PTM':
        # PTM groups don't have comparison suffix
        full_group_name = base_group_name
    else:
        # TP groups have comparison suffix
        full_group_name = f"{base_group_name} ({comparison_label})"
    
    return (full_group_name, position)


def _apply_significance_group_colors(payload: Dict, colors: Dict[str, str], p_cutoff: float, fc_cutoff: float, curtain_type: str):
    """Apply colors to significance group patterns."""
    is_ptm = curtain_type.upper() == 'PTM'
    
    diff_form = payload.get('differentialForm', {})
    comparison_select = diff_form.get('_comparisonSelect', ['1'])
    
    if isinstance(comparison_select, str):
        comparison_select = [comparison_select]
    
    for comparison in comparison_select:
        for position_key, color in colors.items():
            # Build the actual group name with values
            parts = []
            position_parts = position_key.split('FC')
            
            # Handle P-value part
            if 'P-value > ' in position_parts[0]:
                parts.append(f"P-value > {p_cutoff}")
            elif 'P-value <= ' in position_parts[0]:
                parts.append(f"P-value <= {p_cutoff}")
            
            # Handle FC part  
            if len(position_parts) > 1:
                fc_part = position_parts[1].strip()
                if fc_part.startswith('> '):
                    parts.append(f"FC > {fc_cutoff}")
                elif fc_part.startswith('<= '):
                    parts.append(f"FC <= {fc_cutoff}")
            
            base_group_name = ";".join(parts)
            
            # Add comparison suffix based on curtain type
            if is_ptm:
                full_group_name = base_group_name
            else:
                full_group_name = f"{base_group_name} ({comparison})"
            
            payload['settings']['colorMap'][full_group_name] = color


def _search_proteins_frontend_style(payload: Dict, protein_ids: List[str], raw_df: pd.DataFrame = None,
                                   primary_id_col: str = None) -> List[str]:
    """Replicate exact frontend protein ID to Primary ID conversion."""
    result = []
    is_curtain_ptm = payload.get('curtain_type') == 'PTM'
    
    # Process each protein ID individually (frontend processes one at a time)
    for protein_input in protein_ids:
        # Split semicolon-separated values like frontend
        for protein_part in protein_input.split(';'):
            protein_id = protein_part.strip().upper()
            if not protein_id:
                continue
                
            if is_curtain_ptm:
                # CurtainPTM: getPrimaryFromAcc(acc)
                primary_ids = _get_primary_from_acc_ptm(payload, protein_id)
            else:
                # Curtain TP: getPrimaryIDsFromAcc(primaryIDs)
                primary_ids = _get_primary_ids_from_acc_tp(payload, protein_id)
            
            # Add found primary IDs to result
            for primary_id in primary_ids:
                if primary_id not in result:
                    result.append(primary_id)
    
    return result

def _get_primary_from_acc_ptm(payload: Dict, acc: str) -> List[str]:
    """CurtainPTM: getPrimaryFromAcc() exact replication."""
    result = []
    if 'extraData' not in payload or 'data' not in payload['extraData']:
        return result
    
    accession_to_primary_ids = payload['extraData']['data'].get('accessionToPrimaryIDs', {})
    if acc in accession_to_primary_ids:
        primary_list = accession_to_primary_ids[acc]
        if isinstance(primary_list, list):
            for primary_id in primary_list:
                if primary_id not in result:
                    result.append(primary_id)
        elif isinstance(primary_list, str):
            result.append(primary_list)
    
    return result

def _get_primary_ids_from_acc_tp(payload: Dict, primary_ids: str) -> List[str]:
    """Curtain TP: getPrimaryIDsFromAcc() exact replication."""
    result = []
    if 'extraData' not in payload or 'data' not in payload['extraData']:
        return result
    
    primary_ids_map = payload['extraData']['data'].get('primaryIDsmap', {})
    if primary_ids in primary_ids_map:
        primary_dict = primary_ids_map[primary_ids]
        if isinstance(primary_dict, dict):
            for acc in primary_dict.keys():
                if acc not in result:
                    result.append(acc)
    
    return result


def create_search_group_from_regex(payload: Dict, group_name: str, pattern: str,
                                 raw_df: pd.DataFrame, search_column: str,
                                 primary_id_col: str, color: str = None,
                                 specific_comparison_label: str = None) -> Dict:
    """
    Create a search group by matching a regex pattern in a specific column.
    
    Args:
        payload (Dict): The payload to modify.
        group_name (str): Name of the search group.
        pattern (str): Regex pattern to match.
        raw_df (pd.DataFrame): Raw data DataFrame to search.
        search_column (str): Column name to search in.
        primary_id_col (str): Column name containing primary protein IDs.
        color (str, optional): Hex color code for the group.
        specific_comparison_label (str, optional): Specific comparison label override.
    
    Returns:
        Dict: Updated payload with the new search group containing matched proteins.
        
    Example:
        >>> # Find all proteins with "kinase" in their names
        >>> payload = create_search_group_from_regex(
        ...     payload,
        ...     "Kinase proteins",
        ...     r"(?i)kinase",
        ...     raw_df,
        ...     "Protein.names",
        ...     "Protein.IDs"
        ... )
    """
    import re
    
    if isinstance(raw_df, str):
        raw_df = pd.read_csv(raw_df, sep='\t' if raw_df.endswith('.tsv') else ',')
    
    if search_column not in raw_df.columns:
        raise ValueError(f"Search column '{search_column}' not found in DataFrame")
    
    # Find matches using regex
    mask = raw_df[search_column].str.contains(
        pattern, case=False, na=False, regex=True
    )
    matched_proteins = raw_df[mask][primary_id_col].tolist()
    
    if not matched_proteins:
        print(f"Warning: No proteins found matching pattern '{pattern}' in column '{search_column}'")
        return payload
    
    return add_search_group(payload, group_name, matched_proteins, color, specific_comparison_label)


def update_search_group_colors(payload: Dict, color_mapping: Dict[str, str]) -> Dict:
    """
    Update colors for existing search groups.
    
    Args:
        payload (Dict): The payload to modify.
        color_mapping (Dict[str, str]): Mapping of group names to new hex colors.
    
    Returns:
        Dict: Updated payload with new colors assigned.
        
    Example:
        >>> color_mapping = {
        ...     "Kinases": "#FF0000",
        ...     "Phosphatases": "#00FF00"
        ... }
        >>> payload = update_search_group_colors(payload, color_mapping)
    """
    if 'settings' not in payload:
        payload['settings'] = {}
    if 'colorMap' not in payload['settings']:
        payload['settings']['colorMap'] = {}
    
    for group_name, color in color_mapping.items():
        if group_name in payload.get('selectionsName', []):
            payload['settings']['colorMap'][group_name] = color
        else:
            print(f"Warning: Search group '{group_name}' not found in payload")
    
    return payload


def remove_search_group(payload: Dict, group_name: str) -> Dict:
    """
    Remove a search group from the payload.
    
    Args:
        payload (Dict): The payload to modify.
        group_name (str): Name of the search group to remove.
    
    Returns:
        Dict: Updated payload with the search group removed.
    """
    # Remove from operation names
    if 'selectionsName' in payload and group_name in payload['selectionsName']:
        payload['selectionsName'].remove(group_name)
    
    # Remove from color map
    if 'settings' in payload and 'colorMap' in payload['settings']:
        payload['settings']['colorMap'].pop(group_name, None)
    
    # Remove from selectedMap
    if 'selectionsMap' in payload:
        for protein_id in list(payload['selectionsMap'].keys()):
            if group_name in payload['selectionsMap'][protein_id]:
                del payload['selectionsMap'][protein_id][group_name]
                # Remove protein entry if no groups remain
                if not payload['selectionsMap'][protein_id]:
                    del payload['selectionsMap'][protein_id]
    
    return payload


def get_search_group_summary(payload: Dict) -> Dict[str, Dict]:
    """
    Get a summary of all search groups in the payload.
    
    Args:
        payload (Dict): The payload to analyze.
    
    Returns:
        Dict[str, Dict]: Summary with group names as keys and info as values containing:
            - 'protein_count': Number of proteins in the group
            - 'color': Assigned color
            - 'proteins': List of protein IDs (up to 10 shown)
            - 'is_comparison_specific': Whether group is comparison-specific
            
    Example:
        >>> summary = get_search_group_summary(payload)
        >>> print(summary)
        {
            'Kinases': {
                'protein_count': 45,
                'color': '#FF6B6B',
                'proteins': ['P12345', 'Q67890', '...'],
                'is_comparison_specific': False
            }
        }
    """
    summary = {}
    
    if 'selectionsName' not in payload:
        return summary
    
    for group_name in payload['selectionsName']:
        # Count proteins in this group
        protein_count = 0
        proteins = []
        
        if 'selectionsMap' in payload:
            for protein_id, groups in payload['selectionsMap'].items():
                if group_name in groups and groups[group_name]:
                    protein_count += 1
                    if len(proteins) < 10:  # Show up to 10 proteins
                        proteins.append(protein_id)
        
        # Get color
        color = "#808080"  # Default gray
        if ('settings' in payload and 'colorMap' in payload['settings'] 
            and group_name in payload['settings']['colorMap']):
            color = payload['settings']['colorMap'][group_name]
        
        # Check if comparison-specific
        is_comparison_specific = '(' in group_name and group_name.endswith(')')
        
        summary[group_name] = {
            'protein_count': protein_count,
            'color': color,
            'proteins': proteins,
            'is_comparison_specific': is_comparison_specific
        }
    
    return summary
