from pkg import Downloader

with open('api_key.txt', 'r', encoding='utf-8') as file:
    api_key = file.read()

dwn = Downloader(api_key)

dwn.download_to_file('materials.materialarchive')
