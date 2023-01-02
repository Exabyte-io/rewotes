from mp_api.client import MPRester


with MPRester(api_key='<enter your api key>') as mpr:

    # for a single material
    thermo_doc = mpr.thermo.get_data_by_id('mp-1103503')

    # for many materials, it's much faster to use
    # the `search` method, where additional material_ids can 
    # be added to this list
    thermo_docs = mpr.thermo.search(material_ids=['mp-1103503'])

