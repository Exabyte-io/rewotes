from mp_api.client import MPRester

with open('api_key.txt', 'r', encoding='utf-8') as file:
    api_key = file.read()

with MPRester(api_key) as mpr:
    docs = mpr.summary.search(material_ids=['mp-149'])
    print(docs[0])
