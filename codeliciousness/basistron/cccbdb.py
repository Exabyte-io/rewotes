# -*- coding: utf-8 -*-
import requests
from functools import partialmethod
from urllib.parse import urljoin

from typing import Dict, Any, Optional

from bs4 import BeautifulSoup
from bs4.element import Tag

from basistron import utils
from basistron import parser

log = utils.get_logger(__name__)

def _borrowed_headers(referer):
    return {
        'Host': 'cccbdb.nist.gov',
        'Connection': 'keep-alive',
        'Content-Length': '26',
        'Pragma': 'no-cache',
        'Cache-Control': 'no-cache',
        'Origin': 'http://cccbdb.nist.gov',
        'Upgrade-Insecure-Requests': '1',
        'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/58.0.3029.110 Safari/537.36',
        'Content-Type': 'application/x-www-form-urlencoded',
        'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,*/*;q=0.8',
        'Referer': referer,
        'Accept-Encoding': 'gzip, deflate',
        'Accept-Language': 'en-CA,en-GB;q=0.8,en-US;q=0.6,en;q=0.4',
    }

class Cccbdb:
    """Wrapper around interactivity with the CCCBDB and
    filesystem caching of CCCBDB data."""

    TIMEOUT = 30000
    EXT_FMT = "{}x.asp".format
    BASE_URL = "http://cccbdb.nist.gov"
    FORM_EXT = "getform"
    DUMP_EXT = "carttabdump"
    EXPT_EXT = "exp1"

    POL_CALC = "polcalc1"

    def get_calculated_data(self, formula: str, property: str = POL_CALC):
        soup = self.post(
            self.FORM_EXT,
            referer=property,
            data={
                "formula": formula,
                "submit1": "Submit",
            }
        )
        print(soup)
        tables = soup.find_all("table")
        print(len(tables), [len(table) for table in tables])
        import pprint
        # pprint.pprint(tables[0])
        # pprint.pprint(tables[1])
#        tables = soup.find_all("table")
#        log.info(f"found {len(tables)} tables")
#        parse = parser.TableParser()
#        parse.feed(str(soup))
#        for table in parse.tables:
#            log.info(f"table({len(table)},{len(table[0])})")
#            print(table)


    def get_experimental_data(self, formula: str, path: str = EXPT_EXT):
        [form] = self.get(path).find_all("form")
        reduced = self.reduce_form(form)
        data = {inp["name"]: inp["value"] for inp in reduced["inputs"]}
        resp = self.post(reduced["action"], data=data)
        tables = resp.find_all("table")
        log.info(f"found {len(tables)} tables")
        parse = parser.TableParser()
        parse.feed(str(resp))
        print(parse.tables)
    
    @staticmethod
    def reduce_form(form: Tag) -> Dict[str, Any]:
        return {
            "action": form.attrs.get("action").lower(),
            "method": form.attrs.get("method", "get").lower(),
            "inputs": [
                {
                    "type": inp.attrs.get("type", "text"),
                    "name": inp.attrs.get("name"),
                    "value": inp.attrs.get("value"),
                }
                for inp in form.find_all("input")
            ]
        }

    def __init__(self):
        self.session = requests.Session()

    def _make_request(self, method: str, path: str, referer: Optional[str] = None, **kwargs):
        url = urljoin(self.BASE_URL, self.EXT_FMT(path))
        headers = kwargs.get("headers", {})
        if referer is not None:
            headers.update(_borrowed_headers(url))
            kwargs["allow_redirects"] = False
        kwargs["headers"] = headers
        log.info(f"calling {method} {url}")
        res = self.session.request(method, url, timeout=self.TIMEOUT, **kwargs)
        res.raise_for_status()
        log.info(f"status code {res.status_code}")
        if referer is not None and res.status_code == 302:
            url = urljoin(self.BASE_URL, self.EXT_FMT(referer.replace("1", "2")))
            log.info(f"redirect url get {url}")
            res = self.session.request("get", url, timeout=self.TIMEOUT)
            res.raise_for_status()
            log.info(f"status code {res.status_code}")
        return BeautifulSoup(res.content, "html.parser")

    get = partialmethod(_make_request, "get")
    post = partialmethod(_make_request, "post")

        
if __name__ == "__main__":
    c = Cccbdb()
    c.get_calculated_data("CH4")