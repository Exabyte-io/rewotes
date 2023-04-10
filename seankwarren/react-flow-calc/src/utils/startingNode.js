const startingNode = {
    id: 'instructions',
    type: 'default',
    position: { x: 250, y: 150 },
    selected: true,
    data: { label: 
    `1. Drag and drop a node onto the flowchart\n` +
    `2. Connect nodes using the handles on each node. Ensure to connect an outputs (black) to input (white).\n` +
    `3. Edit input values and select functions from the dropdown to see the updated results in the output nodes.\n` +
    `4. View the JSON representation of nodes and edges on the right.` },
}

export default startingNode;