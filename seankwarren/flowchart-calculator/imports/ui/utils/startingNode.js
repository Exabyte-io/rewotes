const startingNode = {
    id: 'instructions',
    type: 'default',
    position: { x: 100, y: 150 },
    selected: true,
    data: { label: 
    `1. Drag and drop nodes onto the flowchart\n` +
    `2. Connect nodes using the handles on each node. Ensure to connect outputs (black) to inputs (white).\n` +
    `3. Edit input values and select functions to see the updated outputs.\n` +
    `4. View the JSON representation of nodes and edges on the right.` },
}

export default startingNode;