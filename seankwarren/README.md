# Flowchart Calculator

Flowchart Calculator is a Meteor/React application that allows users to create flowcharts representing mathematical operations, functions, and logical comparisons. It's built using React, Meteor, and MongoDB and the React Flow library.

## Features

- Drag and drop interface for adding nodes to the flowchart
- Live calculation and display of results
- Support for arithmetic operations (+, -, *, /, ^), functions (sin, cos, tan, exp), and logical comparisons (>, <, >=, <=, ==, !=)
- View real-time syntax-highlighted JSON representation of nodes and edges
- Interactive flowchart editor powered by React Flow
- Support for saving and loading flowcharts
- Dark mode 🎉

## Usage

- Drag and drop a custom node (input, operation, comparison, or output) from the buttons panel onto the flowchart area.
- Connect nodes using the handles on each node. Ensure to connect the output handle of one node to the input handle of another.
- Edit input nodes or operation/comparison nodes to see the updated output values in the output nodes.
- View the JSON representation of nodes and edges in the right pane.
- Name and save flowcharts and reload them

## Setup

Install dependencies using `npm install`.

Run the development server using `meteor`.
> will launch on <a href="http://localhost:3000/">localhost:3000/</a>

Run tests with `npm run cypress:open`


## Front-end Environment Structure:
    imports
    │
    ├── api
    │   └── flows.js                            // Defines collection on MongoDB
    └── ui
        ├── App.jsx                             // Root component for React
        ├── components                          // All React components
        │   ├── customNodes                     // Custom node types for React Flow
        │   │   ├── BinaryNode.jsx
        │   │   ├── ComparisonNode.jsx
        │   │   ├── InputNode.jsx
        │   │   ├── OutputNode.jsx
        │   │   ├── UnaryNode.jsx
        │   │   ├── index.js                    // Allows for destructured import
        │   │   └── nodeTypes.jsx               // Defines node types for React Flow
        │   ├── mainPage                        // Main page elements
        │   │   ├── FlowchartCalculator.jsx     // Main component containing the app's state
        │   │   ├── FlowchartCanvas.jsx         // Contains the React Flow component
        │   │   ├── JSONViewer.jsx              // Contains the JSON panel
        │   │   └── NodeButtons.jsx             // Contains the drag and drop node buttons
        │   └── reusable                        // Generic components that can be reused anywhere
        │       ├── DarkModeContext.jsx
        │       ├── DarkModeSwitch.jsx          
        │       ├── DraggableButton.jsx         
        │       └── ResizablePane.jsx           
        ├── hooks                               // Hooks for maintaining state
        │   ├── useDraggable.jsx                // Handles Drag & Drop functionality
        │   ├── useInputField.jsx               // Handles InputField state change
        │   ├── useLocalFlowData.jsx            // Handles all state for React Flow
        │   ├── usePrevious.jsx                 // Stores previous value of a state
        │   ├── useRemoteFlowData.jsx           // Handles interaction with flows on DB
        │   ├── useResizeable.jsx               // Stores state for SplitPane
        │   ├── useToggleable.jsx               // Handles state change for any toggleable state
        │   └── useUpdateOutputNodes.jsx        // Contains logic for updating outputs on state change
        ├── styles
        │   ├── a11y_light.css                  // Syntax highlighting styling
        │   └── index.css                       // Main page styling
        └── utils
            ├── calculate.js                    // Recursive calculation function
            ├── createNode.js                   // Node creation logic
            ├── getNodeColor.js                 // Access node colors
            ├── operationDef.js                 // Defines operations for the different node type
            ├── performOperation.js             // Performs a single operation
            └── startingNode.js                 // Contains the starting instructions node






## Complete User Stories:
- [x] As an end user, I should be able to perform calculations using these 4 binary operators (+, -, *, /) and 6 comparison operators (<, >, >=, <=, ==, !=).
- [x] As an end user, when I operate on two values, I should be able to read the output as a number in an output node.
- [x] As an end user, when I compare two values, I should be able to read the output as 'true' or 'false' in an output node.
- [x] As an end user, when I operate on or compare two values, the top-most edge should always come first in the operation. e.g $3\over4$ vs $4\over3$
- [x] As an end user, when I add any nodes or edges, the results in all output nodes should update.
- [x] As an end user, the result should update on any change to the input values or selected operators that I make.
- [x] As an end user, I should be able to construct flowcharts of arbitrary depth and number of output cells.
- [x] As an end user, I should be able to drag desired nodes onto the flowchart, and they should be placed at the dropped position.
- [x] As an end user, I should be able to edit the flowchart and have the JSON update accordingly.
- [x] As an end user, I should be able to use common unary operators such as sin, cos, tan, and e<sup>x</sup>.
- [x] As an end user, I should be able to use the y<sup>x</sup> binary operator, that takes two inputs and calculates the exponential.
- [x] As an end user, node types should be quickly and easily distinguished by shape and color.
- [x] As an end user, I should see some nodes on screen when I launch the application that explain how to use the app.
- [x] As an end user, I should be able to easily read the JSON representation with syntax highlighting.
- [x] As an end user, I should be able to click a 'Save Flow' button to save the JSON datastructure to a MongoDB instance.
- [x] As an end user, I should be able to select a saved flow from a dropdown list and load it in the flowchart and JSON viewers.
- [x] As an end user, I should be able to click on a button called 'clear' to clear the flowchart.
- [x] As a development user, I should be able to automatically test all of the binary, unary, and comparison operations using cypress.
- [x] As a development user, I should be able to automatically test all basic user interactions, including toggling darkmode, __resizing the pane__, clearing the flowchart, dragging and dropping nodes, and drawing edges using cypress.
- [x] As a development user, I should be able to automatically test the save and load procedure using cypress.

### TODO User Stories:
- [ ] As a development user, I should be able to automatically test a complex flowchart like the tanh() function using cypress.
- [ ] As a development user, I should be able to automatically test that the rendered JSON is updated on change, and accurate using cypress.
- [ ] As a development user, I should be able to unit test each react component and hook using jest.