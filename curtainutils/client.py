import json
from copy import deepcopy
from typing import List

import pandas as pd
import requests
import seaborn as sns
from matplotlib import pyplot as plt

from curtainutils.common import curtain_base_payload


class CurtainUniprotData:
    def __init__(self, uniprot):
        self.accMap = {i[0]: i[1] for i in uniprot["accMap"]["value"]}
        self.dataMap = {i[0]: i[1] for i in uniprot["dataMap"]["value"]}
        self.db = self._db_to_df(uniprot["db"])
        if "organism" in uniprot:
            self.organism = uniprot["organism"]
        else:
            self.organism = None
        if "geneNameToAcc" in uniprot:
            self.geneNameToAcc = uniprot["geneNameToAcc"]
        else:
            self.geneNameToAcc = {}
        if "geneNameToPrimary" in uniprot:
            self.geneNameToPrimary = uniprot["geneNameToPrimary"]
        else:
            self.geneNameToPrimary = {}
        self.results = uniprot["results"]

    def _db_to_df(self, db: dict):
        data = []
        for i in db["value"]:
            data.append(i[1])
        return pd.DataFrame(data)

    def get_uniprot_data_from_pi(self, primary_id: str) -> pd.Series | None:
        if primary_id in self.accMap:
            acc_match_list = self.accMap[primary_id]
            if acc_match_list:
                if type(acc_match_list) == str:
                    acc_match_list = [acc_match_list]
                for acc in acc_match_list:
                    if acc in self.dataMap:
                        ac = self.dataMap[acc]
                        filter_db = self.db[self.db["From"] == ac]
                        if not filter_db.empty:
                            return filter_db.iloc[0]
        return None

class CurtainClient:
    def __init__(self, base_url: str, api_key: str = ""):
        self.base_url = base_url
        self.request_session = requests.Session()
        self.refresh_token = ""
        self.access_token = ""
        self.api_key = api_key

    def get_data_filter_list(self):
        if self.api_key != "":
            r = requests.get(f"{self.base_url}/data_filter_list/", headers={"X-Api-Key": self.api_key})
        else:
            r = requests.get(f"{self.base_url}/data_filter_list/")
        if r.status_code == 200:
            res = r.json()
            yield res["results"]
            if res['next']:
                while res['next']:
                    if self.api_key != "":
                        r = requests.get(res['next'], headers={"X-Api-Key": self.api_key})
                    else:
                        r = requests.get(res['next'])
                    if r.status_code == 200:
                        res = r.json()
                        yield res["results"]
                    else:
                        raise ValueError(r.text)

        else:
            raise ValueError(r.text)

    def post_curtain_session(self, payload: dict, file: dict, **kwargs) -> str:
        file = {'file': ('curtain-settings.json', json.dumps(file))}
        if self.api_key != "":
            r = requests.post(f"{self.base_url}/curtain/", data=payload, files=file, headers={"X-Api-Key": self.api_key})
        else:
            r = requests.post(f"{self.base_url}/curtain/", data=payload, files=file)
        if r.status_code == 200:
            return r.json()["link_id"]
        else:
            raise ValueError(r.text)

    def curtain_stats_summary(self, last_n_days: int = 30):
        req = requests.get(f"{self.base_url}/stats/summary/{last_n_days}/")
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

        #rotate x axis labels 90 degrees
        for ax in fig.axes:
            plt.sca(ax)
            plt.xticks(rotation=90)

        fig.tight_layout()
        plt.show()

    def download_curtain_session(self, link_id: str):
        link = f"{self.base_url}/curtain/{link_id}/download/token=/"
        req = requests.get(link)
        if req.status_code == 200:
            data = req.json()
            if "url" in data:
                result = requests.get(data["url"])
                if result.status_code == 200:
                    return result.json()
            else:
                return data
        return None

    def download_sessions_list(self, url_list: List[str]):
        for i in url_list:
            result = self.download_curtain_session(i)
            if result:
                with open(f"{i}_de.txt", "wt") as f:
                    f.write(result["processed"])
                with open(f"{i}_raw.txt", "wt") as f:
                    f.write(result["raw"])


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
            sample_cols: List[str], **kwargs) -> dict:
        """Create payload for curtain session"""
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
            primary_id_raw_col:str,
            sample_cols: List[str],
            peptide_col: str,
            acc_col: str,
            position_col: str,
            position_in_peptide_col: str,
            sequence_window_col: str,
            score_col: str,
            **kwargs) -> dict:
        payload = deepcopy(curtain_base_payload)
        with open(de_file, "rt") as f, open(raw_file, "rt") as f2:
            payload["processed"] = f.read()
            payload["raw"] = f2.read()

        payload["differentialForm"]["_foldChange"] = fc_col
        payload["differentialForm"]["_significant"] = p_col
        payload["differentialForm"]["_comparison"] = comp_col
        payload["differentialForm"]["_comparisonSelect"] = comp_select
        payload["differentialForm"]["_primaryIDs"] = primary_id_de_col
        payload["differentialForm"]["_peptideSequence"] = peptide_col
        payload["differentialForm"]["_accession"] = acc_col
        payload["differentialForm"]["_position"] = position_col
        payload["differentialForm"]["_positionPeptide"] = position_in_peptide_col
        payload["differentialForm"]["_score"] = score_col
        payload["differentialForm"]["_sequenceWindow"] = sequence_window_col
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
            primary_id_raw_col:str,
            sample_cols: List[str],
            peptide_col: str = "",
            acc_col: str = "",
            position_col: str = "",
            position_in_peptide_col: str = "",
            sequence_window_col: str = "",
            score_col: str = "",
            **kwargs) -> dict:
        if peptide_col == "" and acc_col == "" and position_col == "" and position_in_peptide_col == "" and sequence_window_col == "" and score_col == "":
            return self.create_curtain_session_payload(de_file, raw_file, fc_col, transform_fc, transform_significant, reverse_fc, p_col, comp_col, comp_select, primary_id_de_col, primary_id_raw_col, sample_cols, **kwargs)
        else:
            return self.create_curtain_ptm_session_payload(de_file, raw_file, fc_col, transform_fc, transform_significant, reverse_fc, p_col, comp_col, comp_select, primary_id_de_col, primary_id_raw_col, sample_cols, peptide_col, acc_col, position_col, position_in_peptide_col, sequence_window_col, score_col, **kwargs)

    def retrieve_curtain_session(self, link_id: str, token: str = ""):
        req = self.request_session.get(f"{self.base_url}/curtain/{link_id}/download/token={token}/")
        if req.status_code == 200:
            res = req.json()
            if "url" in res:
                data_req = self.request_session.get(res["url"])
                if data_req.status_code == 200:
                    return data_req.json()
            else:
                return res
        else:
            raise ValueError(req.text)

def parse_curtain_from_json(json_path: str):
    with open(json_path, "rt") as f:
        data = json.load(f)
    if type(data["settings"]) == str:
        data["settings"] = json.loads(data["settings"])
    print(data["settings"])
    if "version" in data["settings"]:
        if data["settings"]["version"] == "2":
            pass
    else:
        return parse_old_version(data)

    return data

def parse_old_version(json_data: dict):
    if "colormap" not in json_data["settings"]:
        json_data["settings"]["colormap"] = {}
    if "pCutoff" not in json_data["settings"]:
        json_data["settings"]["pCutoff"] = 0.05
    if "log2FCCutoff" not in json_data["settings"]:
        json_data["settings"]["log2FCCutoff"] = 0.6
    if "dataColumns" in json_data["settings"]:
        json_data["settings"]["dataColumns"] = json_data["settings"]["dataColumns"].split(",")

    json_data["settings"]["pCutoff"] = 0.05
    pass

def parse_v2(json_data: dict):
    pass

