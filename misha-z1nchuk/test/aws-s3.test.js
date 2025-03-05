require('dotenv').config();

const fs = require('fs');
const path = require('path');
const AwsUploadFile = require('../lib/CloudStorages/AWS_S3');

const testUploadFolderName = path.join(__dirname, 'fileToUploadTest');



describe('S3 testing', () => {
    beforeAll( () => {
        fs.mkdir(testUploadFolderName, (err) => {
            expect(err).toBe(null);
        });
        fs.appendFile(path.join(testUploadFolderName, 'someFile1.txt'), 'someData', (err) => {
            expect(err).toBe(null);
        });
        fs.appendFile(path.join(testUploadFolderName, 'someFile2.txt'), 'someData', (err) => {
            expect(err).toBe(null);
        });
        fs.appendFile(path.join(testUploadFolderName, 'someFile3.txt'), 'someData', (err) => {
            expect(err).toBe(null);
        });

        jest.setTimeout(10000);
    });


    test('S3 wrong file path', async () => {
        try {
            await AwsUploadFile(
                path.join(testUploadFolderName, 'someFile3.txt'),
                'st',
                process.env.AWS_BUCKET_NAME,
            );
        } catch (e) {
            expect(e).toBe('File not exists');
        }
    });

    test('S3 success upload, if correct file path', async () => {
        try {
            const res = await AwsUploadFile(
                path.join(testUploadFolderName, 'someFile3.txt'),
                path.join(testUploadFolderName, 'someFile3.txt'),
                process.env.AWS_BUCKET_NAME,
            );
            expect(res).toHaveProperty('Location');
            expect(res).toHaveProperty('ETag');
            expect(res).toHaveProperty('Bucket');
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


