from unittest import TestCase
from curtainutils.client import parse_curtain_from_json, CurtainClient


class Test(TestCase):
    def test_parse_curtain_from_json(self):
        json_path = r"C:\Users\Toan Phung\curtain_media\files\curtain_upload\0a0d33de-478b-4654-8cec-811a86400062.json"
        parse_curtain_from_json(json_path)

    def test_get_curtain_filter_list(self):
        data_list = CurtainClient("https://celsus.muttsu.xyz").get_data_filter_list()
        for data in data_list:
            continue

    def test_get_curtain_download(self):
        data = CurtainClient("https://celsus.muttsu.xyz").download_curtain_session("f4b009f3-ac3c-470a-a68b-55fcadf68d0f")
        print(data["settings"])
        print(data["differentialForm"])
        print(data["rawForm"])
        print(data["annotatedData"])
        for e in data["extraData"]:
            print(e)
        for k in data:
            print(k)
        #print(data["processed"])
        #print(data["raw"])