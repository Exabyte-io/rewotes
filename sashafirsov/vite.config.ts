/// <reference types='vitest' />
import { defineConfig } from 'vite';
import react from '@vitejs/plugin-react';
import { nxViteTsPaths } from '@nx/vite/plugins/nx-tsconfig-paths.plugin';
import { configDefaults } from 'vitest/config';

export default defineConfig({
    root: __dirname,
    cacheDir: './node_modules/.vite/.',

    server: {
        port: 4200,
        host: 'localhost'
    },

    preview: {
        port: 4300,
        host: 'localhost'
    },

    plugins: [react(), nxViteTsPaths()],

    // Uncomment this if you are using workers.
    // worker: {
    //  plugins: [ nxViteTsPaths() ],
    // },

    build: {
        outDir: './dist/sashafirsov',
        reportCompressedSize: true,
        commonjsOptions: {
            transformMixedEsModules: true
        }
    },

    test: {
        globals: true,
        cache: {
            dir: './node_modules/.vitest'
        },
        environment: 'jsdom',
        include: ['src/**/*.{test,spec}.{js,mjs,cjs,ts,mts,cts,jsx,tsx}'],

        reporters: ['default'],
        coverage: {
            reportsDirectory: './coverage/sashafirsov',
            provider: 'v8',
            exclude: [...configDefaults.exclude, 'e2e/*', 'src/main.tsx']
        }
    }
});
