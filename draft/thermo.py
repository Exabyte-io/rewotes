from mp_api.client import MPRester
import secret

with MPRester(api_key=secret.API_KEY) as mpr:

    # for a single material
    # thermo_doc = mpr.thermo.get_data_by_id('mp-1103503')      #   DOES NOT WORK (WHY?)

    # for many materials, it's much faster to use
    # the `search` method, where additional material_ids can 
    # be added to this list
    thermo_docs = mpr.thermo.search(material_ids=['mp-1103503'])
    
    print(thermo_docs)
