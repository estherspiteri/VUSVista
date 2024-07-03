import json
from typing import List

import requests
from flask import current_app
from requests import RequestException

from server.responses.internal_response import InternalResponse


def get_consequences_for_new_vus(hgvs: List[str]) -> InternalResponse:
    simplified_hgvs_dict = {}

    for h in hgvs:
        simplified_hgvs = h.split(' ')[0]
        simplified_hgvs_dict[simplified_hgvs] = h

    data = {"hgvs_notations": list(simplified_hgvs_dict.keys())}

    try:
        headers = {"Content-Type": "application/json", "Accept": "application/json"}
        res = requests.post("https://rest.ensembl.org/vep/human/hgvs", headers=headers,
                      data=json.dumps(data))

    except RequestException as e:
        current_app.logger.error(f'Failed to connect to Ensembl Service: {e}')
        # send service unavailable status code
        return InternalResponse(None, 503, e.response)

   # could not parse the HGVS notation
    if res.status_code == 400:
        # TODO: show error msg on front-end to re-enter hgvs
        return InternalResponse({'consequences_dict': {}}, 200)
    if res.status_code != 200:
        current_app.logger.error(f'Ensembl Service failed: {res.reason}')
        return InternalResponse(None, res.status_code, res.reason)
    else:
        res_json = res.json()

        hgvs_dict = {}

        for r in res_json:
            unsimplified_hgvs = simplified_hgvs_dict[r["id"]]
            hgvs_dict[unsimplified_hgvs] = r["most_severe_consequence"]

        return InternalResponse({'consequences_dict': hgvs_dict}, 200)
