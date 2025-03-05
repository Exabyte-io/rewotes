/** @type {import('ts-jest').JestConfigWithTsJest} */
module.exports = {
  preset: 'ts-jest',
  testEnvironment: 'jsdom',
  transform: {
    '^.+\\.jsx?$': 'ts-jest',
  },
  setupFilesAfterEnv: ['./setupTests.ts'],
};
