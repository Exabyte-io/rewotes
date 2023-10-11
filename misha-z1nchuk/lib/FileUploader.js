/* eslint-disable no-console */
/* eslint-disable no-throw-literal */
const path = require('path');
const fs = require('fs');
const readdir = require('recursive-readdir');
const Workers = require('piscina');
const os = require('os');

const rootFolder = path.resolve(__dirname, '../');

// get array of path to files, that we need to upload
function getFiles(dirPath) {
    return fs.existsSync(dirPath) ? readdir(dirPath) : [];
}

class FileUploader {
    workerPool; // worker threads pool, with queue

    constructor({ threads = undefined }) {
        this.workerPool = new Workers({
            filename: path.resolve(__dirname, 'thread.js'), // path to worker instructions
            maxThreads: threads || os.cpus().length, // config max amount of threads
        });
    }

    async uploadContent(uploadFolder) {
        if (fs.existsSync(path.resolve(__dirname, uploadFolder))) {
            const filesToUpload = await getFiles(path.resolve(__dirname, uploadFolder));

            await Promise.all(filesToUpload.map(async (file) => {
                try {
                    const Key = file.replace(`${rootFolder}/`, '');

                    console.log(`uploading: [${Key}]`);
                    await this.workerPool.run({
                        file,
                        Key,
                    });

                    console.log(`remaining files amount :${this.workerPool.queueSize}`);
                } catch (e) {
                    console.log(e);
                    throw e;
                }
            }));
        } else {
            throw 'Folder does not exist';
        }
    }
}
module.exports = {
    FileUploader,
    getFiles,
};
