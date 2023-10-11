/* eslint-disable no-console */
/* eslint-disable consistent-return */
require('dotenv').config();
const AwsUploadFile = require('./CloudStorages/AWS_S3');

module.exports = async ({ Key, file }) => {
    try {
        await AwsUploadFile(Key, file, process.env.AWS_BUCKET_NAME);

        return 'done';
    } catch (e) {
        console.log(e);
    }
};
