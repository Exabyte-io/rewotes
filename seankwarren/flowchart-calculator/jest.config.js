module.exports = {
    testEnvironment: 'jsdom',
    setupFilesAfterEnv: ['<rootDir>/imports/setupTests.js'],
    testPathIgnorePatterns: ["/node_modules/", "/tests/"],
    moduleNameMapper: {
        '^meteor/(.*)': '<rootDir>/__mocks__/meteor/$1.js',
        // "\\.(css|less|scss|sass)$": '<rootDir>/__mocks__/styleMock.js',
    },
    transform: {
        "^.+\\.jsx?$": "babel-jest",
    }
};