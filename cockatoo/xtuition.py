import os
import requests
import cockatoo

base_uri = 'http://xtuition.org/api'

def _auth():
    if 'XTUITION_TOKEN' not in os.environ:
        raise(ValueError('Please set XTUITION_TOKEN environment variable'))
    return {'Authorization': 'Bearer {}'.format(os.environ['XTUITION_TOKEN'])}

def _fetch_json(endpoint, payload=None):
    r = requests.get(endpoint, headers=_auth(), params=payload)
    return r

def fetch_screen(id):
    endpoint = base_uri + '/screen/' + str(id) + '/cockatoo'
    r = _fetch_json(endpoint)
    return cockatoo.screen.loads(r.text)

def fetch_cocktail(id):
    endpoint = base_uri + '/cocktail/' + str(id) + '/cockatoo'
    r = _fetch_json(endpoint)
    return cockatoo.screen._parse_cocktail_json(r.json())
