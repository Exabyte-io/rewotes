const { program } = require('commander');

program
    .name('Multi-threaded-s3-file-uploader')
    .description('CLI')
    .version('0.0.1');

program
    .option('--threads <int>', 'max amount of threads to use')
    .option('-f, --folder <string>', "path to folder for uploading, from lib folder");

program.parse();

const options = program.opts();
const threads = options.threads || undefined;
const pathFolder = options.folder || undefined;

module.exports = {
    program,
    threads: parseInt(threads, 10),
    pathFolder,
};
