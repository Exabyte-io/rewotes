# -*- coding: utf-8 -*-
"""
A CCCBDB scraping client. Borrows much inspiration
from https://github.com/marcelo-mason/cccbdb-calculation-parser
but the repo is dormant and non-functional.
"""
import sys
from functools import partialmethod
from typing import Any, Dict, List, Tuple
from urllib.parse import urljoin

import pandas as pd
import requests
from bs4 import BeautifulSoup
from bs4.element import Tag
from requests import HTTPError, ReadTimeout

from basistron import parser, utils

log = utils.get_logger(__name__)

def _inspected_headers(referer: str) -> Dict[str, str]:
    headers = {
        'Accept': (
            'text/html,application/xhtml+xml,application/xml;'
            'q=0.9,image/avif,image/webp,image/apng,*/*;'
            'q=0.8,application/signed-exchange;v=b3;q=0.9'
        ),
        'Accept-Encoding': 'gzip, deflate, br',
        'Accept-Language': 'en-US,en;q=0.9',
        'Cache-Control': 'max-age=0',
        'Connection': 'keep-alive',
        'Content-Length': '26',
        'Content-Type': 'application/x-www-form-urlencoded',
        'Host': 'cccbdb.nist.gov',
        'Origin': 'https://cccbdb.nist.gov',
        'sec-ch-ua': (
            '"Chromium";v="94", "Google Chrome";'
            'v="94", ";Not A Brand";v="99"'
        ),
        'sec-ch-ua-mobile': '?0',
        'sec-ch-ua-platform': 'Windows',
        'Sec-Fetch-Dest': 'document',
        'Sec-Fetch-Mode': 'navigate',
        'Sec-Fetch-Site': 'same-origin',
        'Sec-Fetch-User': '?1',
        'Upgrade-Insecure-Requests': '1',
        'User-Agent': (
            'Mozilla/5.0 (Windows NT 10.0; Win64; x64) '
            'AppleWebKit/537.36 (KHTML, like Gecko) '
            'Chrome/94.0.4606.81 Safari/537.36'
        )
    }
    if referer is not None:
        headers['Referer'] = referer
    return headers

class Cccbdb:
    """Wrapper around interactivity with the CCCBDB and
    filesystem caching of CCCBDB data."""

    TIMEOUT = 10
    BASE_URL = "https://cccbdb.nist.gov"
    FORM_EXT = "getformx.asp"
    DUMP_EXT = "carttabdumpx.asp"

    # experimental properties
    EXPT_POL = "pollistx.asp"
    EXPT_VIB = "expvibs1x.asp"
    EXPT_IE = "xp1x.asp?prop=8"
    # general purpose
    # EXPT_EXT = "exp1x.asp"

    # calculated properties
    POL_CALC = "polcalc1x.asp"
    VIB_FREQ = "vibs1x.asp"
    HOMO_LUMO = "gap1x.asp"
    # extend this for more properties
    # ENERGY = "energy1x.asp"
    # MULLIKEN = "mulliken1x.asp"

    @staticmethod
    def retry_loop(func, *args, **kwargs):
        tries = 0
        while True:
            try:
                tries += 1
                if not tries or not tries % 5:
                    log.info(f"calling {func.__name__} try #{tries}")
                res = func(*args, **kwargs)
                break
            except (HTTPError, ReadTimeout):
                continue
            except KeyboardInterrupt:
                sys.exit()
        return res

    def get_data(self, formula: str, property: str = POL_CALC):
        """Get all the data for a chemical formula and kind of property."""
        res, form = self.get_form(property)
        form_data = self.reduce_form(form, formula)
        headers = _inspected_headers(res.url)
        res = self.submit_form(form_data, headers)
        return self.soup(res)

    def get_form(self, property: str) -> Tuple[requests.Response, Dict[str, Any]]:
        """Get the structure of the form directly from the website."""
        res = self.retry_loop(self.get, property)
        soup = self.soup(res)
        # assert only a single form on the page
        [form] = soup.find_all("form")
        return res, form

    def submit_form(self, form_data: Dict[str, Any], headers: Dict[str, Any]):
        """Submit the form following redirect semantics of the CCCBDB website."""
        log.info('submitting form %s', form_data["data"])
        # post the form data without redirect
        self.retry_loop(
            getattr(self, form_data["method"]),
            form_data["action"],
            data=form_data["data"],
            allow_redirects=False,
            headers=headers,
        )
        # get the data from the redirected URL
        url = headers["Referer"].split("?")[0].replace("1", "2")
        return self.retry_loop(self.get, url)

    @staticmethod
    def reduce_form(form: Tag, formula: str) -> Dict[str, Any]:
        reduced = {
            "action": form.attrs.get("action").lower(),
            "method": form.attrs.get("method", "get").lower(),
            "data": {
                inp.attrs.get("name"): inp.attrs.get("value")
                for inp in form.find_all("input")
            }
        }
        reduced["data"]["formula"] = formula
        return reduced

    @staticmethod
    def soup(response: requests.Response) -> BeautifulSoup:
        return BeautifulSoup(response.content, "html.parser")

    def get_dataframes(self, formula: str, property: str) -> List[pd.DataFrame]:
        soup = self.get_data(formula, property)
        tables = soup.find_all("table")
        log.info(f"found {len(tables)} tables on page")
        dfs = []
        for html in tables:
            parsed = parser.TableParser()
            parsed.feed(str(html))
            df = parsed.to_df()
            if df is not None:
                log.info(f"adding dataframe '{df.name}' of shape {df.shape}")
                dfs.append(df)
            else:
                log.info("skipping table")
        return dfs

    def __init__(self):
        self.session = requests.Session()

    def _make_request(self, method: str, path: str, **kwargs):
        url = urljoin(self.BASE_URL, path)
        log.info(f"calling {method} {url}")
        res = self.session.request(
            method, url, timeout=self.TIMEOUT, **kwargs
        )
        res.raise_for_status()
        log.info(f"status code {res.status_code}")
        return res

    get = partialmethod(_make_request, "get")
    post = partialmethod(_make_request, "post")