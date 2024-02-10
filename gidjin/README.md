# End-to-end Tests (DevOps)

> Ideal candidate: skilled software engineer versed in application infrastructure and DevOps

# Overview

The aim of this task is to create a simple application package (either python or javascript) that includes
complete application testing infrastructure as well as a complete CICD solution using Github workflows.

# Requirements

1. A non-trivial application, e.g. a Flask server with a UI or a React app with a UI with testable components
2. An appropriate end-to-end testing framework implementation (e.g. Cypress) for the application
3. An automated workflow using Github actions to verify that the tests pass

# Expectations

- The application may be relatively simple, this is focused more on application infrastructure and DevOps, but the tests must actually verify functionality
- Correctly passes the tests in automation and displays coverage metrics
- Clean workflow logic

# Timeline

We leave exact timing to the candidate. Must fit Within 5 days total.

# User story

As a developer of this application I can:

- view important coverage metrics of my application
- be aware of the number of tests running/passing when developing

# Notes

- Commit early and often

## React + TypeScript + Vite

This template provides a minimal setup to get React working in Vite with HMR and some ESLint rules.

Currently, two official plugins are available:

- [@vitejs/plugin-react](https://github.com/vitejs/vite-plugin-react/blob/main/packages/plugin-react/README.md) uses [Babel](https://babeljs.io/) for Fast Refresh
- [@vitejs/plugin-react-swc](https://github.com/vitejs/vite-plugin-react-swc) uses [SWC](https://swc.rs/) for Fast Refresh

## Expanding the ESLint configuration

If you are developing a production application, we recommend updating the configuration to enable type aware lint rules:

- Configure the top-level `parserOptions` property like this:

```js
export default {
  // other rules...
  parserOptions: {
    ecmaVersion: 'latest',
    sourceType: 'module',
    project: ['./tsconfig.json', './tsconfig.node.json'],
    tsconfigRootDir: __dirname,
  },
}
```

- Replace `plugin:@typescript-eslint/recommended` to `plugin:@typescript-eslint/recommended-type-checked` or `plugin:@typescript-eslint/strict-type-checked`
- Optionally add `plugin:@typescript-eslint/stylistic-type-checked`
- Install [eslint-plugin-react](https://github.com/jsx-eslint/eslint-plugin-react) and add `plugin:react/recommended` & `plugin:react/jsx-runtime` to the `extends` list

# ToDos

- [x] Create a subfolder according to my GitHub username
- [x] Re-use the content from `../` for `README.md` and modify as needed.
- [x] Contribute code in this folder
- [ ] Put a pull request when done
