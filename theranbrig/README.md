# Materials Designer (Frontend/UX)

> Ideal candidate: skilled front-end developer with UI/UX chops.

# Overview

Create a skeleton IDE (integrated development environment) for materials design. Close to Adobe Dreamweaver (or any other IDE) - when you can change html-markup and simultaneously see the result in another tab. We edit material in XYZ format and view result in 3D.

Front-end developers: use Meteor and React.js and minimalistic UX/UI.

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

# AUTHOR NOTES

For this project I tried to remake a basic version of the Materials Designer using React Three Fiber. This is a React wrapper for Three.js.  It allows for cleaner components as we are allowed to utilize React hooks more easily for state updates.  This basic version consists of an editor and a canvas sections.  When you add the coordinates to the grid they appear instantly.  To add line connections you can simply enter in the ID of the two points and it will map a connecting line between those two points.  The editor right now will save a material when you enter in a name and save it the db.  Once a material is saved, you can select saved materials and edit them just as you would with a new material.  This is a basic working version, but there are more features that I would add given more time.

This app is not perfect as it is just a prototype, but the direction of it is headed in the right direction.  Below are some issues that I think could be addressed as the product is developed out more.

## State Management

Due to time constraints for this version I did not setup Redux. However, the complexity of the growing shared state between components does mean that a state management tool is needed.  In order to pass state to components and make sure that that actions are dispatched properly I would add a store to this.  As of now it would contain materials, new connections, new nodes, and other data that is necessary to make a proper editor and stop the large amount of prop drilling that is currently in the app.  Setting up state would be my first priority.

## Testing

Due to time constraints, and limited prior Meteor use I focused on getting the client to work properly.  However, it should have a full test suite to go along with this.  I would add Unit tests for all of the components to make sure that rendering is done properly.  Also adding tests to the Redux actions and dispatches to make sure they are covered is key.

## File Upload

This was not implemented due to time constraints.  I did research on packages that would allow this, but due to time it would not have been feasible.

## Features To Add

- **Canvas Interactions** - Right now all changes to the model are done through the text editor.  Adding features like making line connections by clicking on two noes would allow for greater UX.  Right now the only canvas feature is that when you hover, it shows which node is selected in the table, so that when you create connections or want to delete them, you know which node is the correct one.

- **Robust Editor** - Features like adding colors to different nodes, easier connection creations, and auto save would be features I would like to implement.  These things would all be easier with more robust state management and more time.

- **User Flow** - The flow should have a screen where you either select a saved material or create a new one, rather than all being on one screen.  This would allow for greater separation of concerns and a better UI.  Right now the editor is a bit cluttered.

- **Form Safety** - Right now there are just some basics for the form safety.  However, making sure that all inputs are safe is a key to a good product.  This of course should be coupled with the backend to make sure that the whole application is safe.
