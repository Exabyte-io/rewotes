import json

import gevent
from gevent import monkey  # pylint: disable=import-error
monkey.patch_all() # patches blocking http I/O by non-blocking variant

import libcloud
from libcloud.storage.types import Provider
from libcloud.storage.providers import get_driver

import Config

class BaseUploader:
    def __init__(self, provider):
        self.provider = provider
        self.config = Config.config
        self.factory = get_driver(provider)
        self.driver = self.initDriver()
        self.container = self.initContainer()

    def initDriver(self):
        return self.factory(
          key=self.config.get('KEY'),
          secret=self.config.get('SECRET'),
          signature_version='4'
        )
    
    def uploadSimple(self, driver, src, dst, meta_data={}):
        extra = {'meta_data': meta_data} # {'acl': 'public-read'}
        obj = driver.upload_object(src, self.container, dst, extra=extra)
        return obj

    def uploadMultiPart(self, driver, src, dst, meta_data={}):
        print('using uploadMultiPart for:', src)
        extra = {'meta_data': meta_data} # {'acl': 'public-read'}
        with open(src, 'rb') as iterator:
            obj = driver.upload_object_via_stream(
                iterator=iterator,
                container=self.container,
                object_name=dst,
                extra=extra
            )
            return obj

    def initContainer(self):
        name = self.config.get('CONTAINER_NAME')
        containers = self.driver.list_containers()
        container = next((el for el in containers if el.name==name), None)
        if not container:
            print('creating container...')
            container = self.driver.create_container(name)
        return container

    def upload(self, src, dst, meta_data={}, filesize=0):
        # use multipart upload for +X MB files
        if (filesize/1024/1024 >= self.config.get('MULTIPART_SIZE_MB')): 
            return self.uploadMultiPart(self.driver, src, dst, meta_data)
        else:
            return self.uploadSimple(self.driver, src, dst, meta_data)

    def saveMeta(self):
        containers = self.driver.iterate_container_objects(self.container)
        meta = {}
        for cont in containers:
            meta[cont.name] = {
                'size': cont.size,
                'last_modified': cont.extra['last_modified'],
                'acl': None,
            }

        with open('meta.json', 'w') as outfile:
            json.dump(meta, outfile, indent=4)
