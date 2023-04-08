import performOperation from './performOperation';

// Function to calculate the value of a node
const calculate = (nodes, edges, sourceHandleId) => {
    // Find the node with the specified ID
    const node = nodes.find((n) => n.id === sourceHandleId.split("-")[0]);
    if (!node) return 0;

    // Calculate value based on node type
    if (node.type === 'inputNode') {
        return node.data.value;
    } else if (node.type === 'binaryNode' || node.type === 'unaryNode' || node.type === 'comparisonNode') {
        // Get all incoming edges for the node
        const inputEdges = edges.filter((e) => e.target === node.id);

        // Sort the input edges based on the position of their target handle
        inputEdges.sort((edge1, edge2) => {
            const edge1Position = edge1.targetHandle.split('-')[2];
            const edge2Position = edge2.targetHandle.split('-')[2];
          
            return edge1Position === 'top' ? -1 : 1;
          });

        // Calculate the values of all incoming nodes
        const inputValues = inputEdges.map((e) => {
            return calculate(nodes, edges, e.sourceHandle);
        });


        // Perform the specified operation on the input values
        const result = performOperation(
            node.data.value, 
            inputValues[0], 
            inputValues[1]
        );
        // console.log(result);
        return result;
    }

    return 0;
}

export default calculate;
