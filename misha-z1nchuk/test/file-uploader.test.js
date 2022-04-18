require('dotenv').config();

const path = require('path');
const { getFiles, FileUploader } = require('../lib/FileUploader');
const fs = require("fs");

const testUploadFolderName = path.join(__dirname, 'fileToUploadTestTwo');



describe('FileUploader module', () => {
    beforeAll(() => {
        fs.mkdir(testUploadFolderName, (err) => {
            expect(err).toBe(null);
        });
        fs.appendFile(path.join(testUploadFolderName, 'someFile1Two.txt'), 'someData', (err) => {
            expect(err).toBe(null);
        });
        fs.appendFile(path.join(testUploadFolderName, 'someFile2Two.txt'), 'someData', (err) => {
            expect(err).toBe(null);
        });
        fs.appendFile(path.join(testUploadFolderName, 'someFile3Two.txt'), 'someData', (err) => {
            expect(err).toBe(null);
        });

        jest.setTimeout(15000);
    });


    test('Get all files path, that are in upload folder', async () => {
        const resultGetFiles = await getFiles(testUploadFolderName);
        const amountOfCreatedFiles = 3;
        expect(resultGetFiles).toBeInstanceOf(Array);
        expect(resultGetFiles.length).toBe(amountOfCreatedFiles);
    });

    test('File Uploader wrong path to folder', async () => {
        try {
            const fileUploader = new FileUploader({});
            await fileUploader.uploadContent('someWrongPath');
        } catch (e) {
            expect(e).toBe('Folder does not exist');
        }
    });

    test('File Uploader correct path to folder', async () => {
        try {
            const fileUploader = new FileUploader({});
            await fileUploader.uploadContent(testUploadFolderName);
        } catch (e) {
            expect(e).toBe(null);
        }
    });


    afterAll(async () => {
        await fs.rmdir(testUploadFolderName, {recursive: true}, (err) => {
            expect(err).toBe(null);
        });
    });
});


