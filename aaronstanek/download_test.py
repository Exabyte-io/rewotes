from pkg import Downloader

with open('api_key.txt', 'r', encoding='utf-8') as file:
    api_key = file.read()

dwn = Downloader('materials.materialarchive', api_key)

elapsed = dwn.download()

print('Time Elapsed:', elapsed)
