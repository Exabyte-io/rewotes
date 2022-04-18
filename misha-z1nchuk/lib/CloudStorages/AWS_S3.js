/* eslint-disable no-throw-literal */
const S3 = require('aws-sdk/clients/s3');
const stream = require('stream');
const fs = require('fs');


if (!process.env.AWS_BUCKET_REGION ||
    !process.env.AWS_ACCESS_KEY ||
    !process.env.AWS_ACCESS_SECRET_KEY) throw "All credentials for S3 is required"


const s3 = new S3({
    region: process.env.AWS_BUCKET_REGION,
    accessKeyId: process.env.AWS_ACCESS_KEY,
    secretAccessKey: process.env.AWS_ACCESS_SECRET_KEY,
});

async function AwsUploadFile(Key, file, Bucket) {
    if (!fs.existsSync(file)) {
        throw 'File not exists';
    }
    const pass = new stream.PassThrough();
    fs.createReadStream(file).pipe(pass);

    return s3.upload({
        Key,
        Bucket,
        Body: pass,
    }).promise();
}
module.exports = AwsUploadFile;
