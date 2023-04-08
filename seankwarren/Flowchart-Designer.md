# Flowchart Calculator

Flowchart Calculator is a React application that allows users to create flowcharts representing mathematical operations and logical comparisons. Users can visually create and connect nodes representing input values, operations, comparisons, and output results.
### Features

- Drag and drop interface for adding nodes to the flowchart
- Custom node types: input, operation, comparison, and output
- Visual connections between nodes using edges
- Automatic calculation and display of output node results based on connected input and operation nodes
- Support for arithmetic operations (+, -, *, /) and logical comparisons (>, <, >=, <=, ==, !=)
- Real-time update of output nodes when input values or operations are changed
- Interactive flowchart editor powered by React Flow

### Use Case

This application can be used as a visual calculator for solving complex mathematical expressions and logical comparisons. Users can create a flowchart by connecting input nodes to operation and comparison nodes, which can then be connected to output nodes to display the results. By adjusting the input values and operations, users can easily explore different scenarios and observe the corresponding output changes in real time.

### How to Run:

install dependecies:
> `npm install`

then launch:
> `npm run dev`
>
> will launch on <a href="http://localhost:5173/">localhost:5173/</a>


Complete User Stories:
- [x] As an end user, I should be able to perform calculations using these 4 binary operators (+, -, *, /) and 6 comparison operators (<, >, >=, <=, ==, !=)
- [x] As an end user, when I operate on two values, I should be able to read the output as a number in an output node
- [x] As an end user, when I compare two values, I should be able to read the output as 'true' or 'false' in an output node
- [x] As an end user, when I operate on or compare two values, the top-most edge should always come first in the operation. e.g $3\over4$ vs $4\over3$
- [x] As an end user, when I add any nodes or edges, the results in all output nodes should update.
- [x] As an end user, the result should update on any change to the input values or selected operators that I make.
- [x] As an end user, I should be able to construct flowcharts of arbitrary depth and number of output cells.
- [x] As an end user, I should be able to drag desired nodes onto the flowchart, and they should be placed at the dropped position.


Potential Improvement User Stories:
- [ ] As an end user, node types should be quickly and easily distinguished by shape and color.
- [ ] As an end user, I should be able to edit the flowchart and have the JSON update accordingly.
- [ ] As an end user, I should be able to edit the JSON and have the flowchart update accordingly.
- [ ] As an end user, I should be able to click a 'Save Flow' button to save the JSON datastructure to a MongoDB instance.
- [ ] As an end user, I should be able to select a saved flow from a dropdown list and load it in the flowchart and JSON viewers.
- [ ] As an end user, I should be able to use the comparison nodes as a gate to decide whether another node should be executed.
- [ ] As an end user, I should be able to use common unary oeprators such as sin, cos, tan, and e<sup>x</sup>.
- [ ] As an end user, I should be able to use the y<sup>x</sup> binary operator, that takes two inputs and calculates the exponential.
- [ ] As an end user, I should be able to attach a GradientNode to any node, and get the gradient of that value with respect to the final result
- [ ] As a user, the input nodes should be exapandable to include multiple values each with their own handle

