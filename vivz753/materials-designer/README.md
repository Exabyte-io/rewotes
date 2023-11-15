# About

## Introduction

This is my skeleton IDE submission for the Materials Designer test as part of the interview process for Mat3ra.

I put my artistic and engineering skills to work in order to forge this masterpiece. It is missing some implementations that would make it ideal (see TODO list), but it has the main functionalities:

- Source Editor
- Focusable menu/sections
- 3D Visualizer
- Importer for .xyz or .txt files
- Interactive inputs for segment lengths

I know the timeline is limited to 5 days, but the last few commits pertain to minor cosmetic changes.

Thanks for checking it out!

## To view

Deployed here [https://materials-designer.vercel.app/](https://materials-designer.vercel.app/)

or

```bash
npm run dev
# or
yarn dev
# or
pnpm dev
```

[http://localhost:3000](http://localhost:3000)

# TODO

- [x] import .xyz & .txt files
- [ ] import .poscar files
- [x] render custom cube geometry
- [x] render cube on explorer
- [x] add angle + segment inputs
- [x] mobile view
- [x] add response to abc length inputs
- [ ] add response to angle inputs
- [ ] add geometries for diff lattice structures
- [x] polish colors, padding, font sizes
- [ ] add modal for adding material
- [x] submit [job application](https://wellfound.com/jobs/1807964-lead-frontend-engineer-react)
- [ ] add info dialog when hovering over sphere
- [ ] write tests

# Materials Designer (Frontend/UX)

> Ideal candidate: skilled front-end developer with UI/UX chops.

# Overview

Create a skeleton IDE (integrated development environment) for materials design, where one can change text and simultaneously see the visual result in another tab. We edit material in XYZ format and view result in 3D.

Front-end developers: use React.js and minimalistic UX/UI.

Pure UI/UX designers: create high fidelity mockups. 

# Requirements

1. Build a general layout with focus areas (eg. toolbar, structure viewer, source viewer)
2. Implement two edit modes:
   - source editor (to edit material in textual representation)
   - visual editor (to adjust material visual representation)
3. Support import  from a file format, (eg. POSCAR, XYZ)

# Expectations

1. up and running application OR high fidelity clickable mockups
2. general IDE layout (e.g. menu, toolboxes)
3. reactive material editor (edit material file → immediately see results on 3D representation)
4. clean and documented code
5. tests

# Timeline

We leave exact timing to the candidate. Must fit Within 5 days total.

# Notes

- use a designated github repo for version control and submission

# Examples

See [Materials Designer](https://github.com/Exabyte-io/materials-designer) repository also.