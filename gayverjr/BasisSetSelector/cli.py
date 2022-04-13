import argparse
import json
parser = argparse.ArgumentParser(description='Parse a JSON file.')
parser.add_argument('json_file', metavar='JSON File', type=str,help='JSON file specifying job information.')

def read_json(json_file):
    with open(json_file, "r") as f:
        data = json.load(f)
    sections = data['Basis Set Selector']
    for section in sections:
        if section['title'] == 'molecule':
            for mol in section['geometry']:
                pass
        elif section['title'] == 'basis':
            pass
        elif section['title'] == 'reference':
            pass
        else:
            raise RunTimeError("Not a valid section.")
    print("Hello world.")

def main():
    args = parser.parse_args()
    read_json(args.json_file)