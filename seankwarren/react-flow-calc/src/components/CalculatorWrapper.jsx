import React, { useState } from 'react';
import FlowchartViewer from './FlowchartViewer';
import JSONViewer from './JSONViewer';
import NodeButtons from './NodeButtons';
import SplitPane, { Pane } from 'split-pane-react';
import 'split-pane-react/esm/themes/default.css';
import 'reactflow/dist/style.css';

const CalculatorWrapper = () => {
    // Set up state for nodes and edges
    const [nodes, setNodes] = useState([]);
    const [edges, setEdges] = useState([]);

    const [draggedNodeType, setDraggedNodeType] = useState(null);

    const [sizes, setSizes] = useState(['60%', '50%']);

    const handleSizeChange = (sizes) => {
        setSizes(sizes);
    };

    const handleDragStart = (e, nodeType) => {
        setDraggedNodeType(nodeType);
        e.dataTransfer.setData('text/plain', nodeType);
    };

    // Function to add a new node
    const addNode = (type, position) => {
        // Generate a unique ID for the node
        const id = `${type}${nodes.length + 1}`;

        // Create a new node object with an ID, type, initial value, and onChange function
        const newNode = {
            id,
            type: `${type}Node`,
            data: {
                value:
                    type === 'input'
                        ? 0
                        : type === 'operation'
                        ? 'add'
                        : type === 'comparison'
                        ? 'greater'
                        : null,
                onChange: (value) => {
                    // Update the value of the node when its input changes
                    setNodes((ns) => {
                        return (ns.map((n) => 
                            n.id === id
                                ? { ...n, data: { ...n.data, value } }
                                : n
                        ))
                    });
                },
            },
            position: position, // Set the position here
        };

        // Add the new node to the list of nodes
        setNodes((ns) => [...ns, newNode]);
    };

    return (
        <SplitPane split="vertical" sizes={sizes} onChange={handleSizeChange}>
            <div className="flowchart-container" style={{ height: '100vh', width: '100%' }}>
                <NodeButtons addNode={addNode} handleDragStart={handleDragStart} />
                <FlowchartViewer
                    nodes={nodes}
                    edges={edges}
                    setNodes={setNodes}
                    setEdges={setEdges}
                    addNode={addNode}
                    draggedNodeType={draggedNodeType}
                    setDraggedNodeType={setDraggedNodeType}
                />
            </div>
            <JSONViewer nodes={nodes} edges={edges} />
        </SplitPane>
    );
};

export default CalculatorWrapper;
