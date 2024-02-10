import { defineConfig } from 'vite'
import react from '@vitejs/plugin-react'

// https://vitejs.dev/config/
export default defineConfig({
  plugins: [react()],
  test: {
    globals: true,
    environment: 'jsdom',
    include: ['./src/**/*.{test,spec}.{ts,tsx}'],
    setupFiles: './tests.setup.ts',
    coverage: {
      provider: 'v8',
      exclude: ['**/*.config.*', '.eslintrc.cjs', '**/main.tsx', '**/vite-env.d.ts'],
      reporter: ['text', 'json-summary', 'json', 'lcov'],
      // Example threshold
      threshold: {
        statements: 90,
        branches: 90,
        functions: 90,
        lines: 90,
      },
    }
  },
})
