/* eslint-disable no-console */
require('dotenv').config();
const { threads, pathFolder } = require('./commander');
const { FileUploader } = require('./FileUploader');

const uploadFolder = pathFolder || '../folder_to_upload';

const fileUploader = new FileUploader({ threads });

fileUploader.uploadContent(uploadFolder).then(() => {
    console.log('task complete');
    process.exit(0);
}).catch((err) => {
    console.error(err);
    process.exit(1);
});
