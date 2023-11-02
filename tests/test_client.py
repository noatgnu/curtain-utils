from unittest import TestCase
from curtainutils.client import parse_curtain_from_json

class Test(TestCase):
    def test_parse_curtain_from_json(self):
        json_path = r"C:\Users\Toan Phung\curtain_media\files\curtain_upload\0a0d33de-478b-4654-8cec-811a86400062.json"
        parse_curtain_from_json(json_path)