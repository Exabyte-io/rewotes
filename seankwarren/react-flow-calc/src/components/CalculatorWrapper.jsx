import React, { useState } from 'react';
import FlowchartViewer from './FlowchartViewer';
import JSONViewer from './JSONViewer';
import NodeButtons from './NodeButtons';
import SplitPane, { Pane } from 'split-pane-react';
import Switch from 'react-switch';
import 'split-pane-react/esm/themes/default.css';
import 'reactflow/dist/style.css';

const CalculatorWrapper = () => {
    // Set up state for nodes and edges
    const [nodes, setNodes] = useState([]);
    const [edges, setEdges] = useState([]);

    // Drag and drop state
    const [draggedNodeType, setDraggedNodeType] = useState(null);

    // Style states
    const [sizes, setSizes] = useState(['60%', '50%']);
    const [isDarkMode, setIsDarkMode] = useState(true);

    const toggleDarkMode = () => {
        setIsDarkMode(!isDarkMode);
    };

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
                        : type === 'binary'
                        ? 'add'
                        : type === 'comparison'
                        ? 'greater'
                        : null,
                onChange: (value) => {
                    // Update the value of the node when its input changes
                    setNodes((ns) => {
                        return ns.map((n) =>
                            n.id === id
                                ? { ...n, data: { ...n.data, value } }
                                : n
                        );
                    });
                },
            },
            position: position, // Set the position here
        };

        // Add the new node to the list of nodes
        setNodes((ns) => [...ns, newNode]);
    };

    return (
        <SplitPane split='vertical' sizes={sizes} onChange={handleSizeChange}>
            <div
                className='flowchart-container'
                style={{ height: '100vh', width: '100%' }}
            >
                <NodeButtons
                    addNode={addNode}
                    handleDragStart={handleDragStart}
                    isDarkMode={isDarkMode}
                />
                <FlowchartViewer
                    nodes={nodes}
                    edges={edges}
                    setNodes={setNodes}
                    setEdges={setEdges}
                    addNode={addNode}
                    draggedNodeType={draggedNodeType}
                    setDraggedNodeType={setDraggedNodeType}
                    isDarkMode={isDarkMode}
                />
            </div>
            <JSONViewer nodes={nodes} edges={edges} isDarkMode={isDarkMode}>
                <Switch
                    checked={isDarkMode}
                    onChange={toggleDarkMode}
                    offColor='#f0f0f0'
                    onColor='#606060'
                    offHandleColor='#ccc'
                    onHandleColor='#444'
                    handleDiameter={30}
                    uncheckedIcon={false}
                    checkedIcon={false}
                    boxShadow='0px 1px 5px rgba(0, 0, 0, 0.6)'
                    // activeBoxShadow="0px 0px 1px 10px rgba(0, 0, 0, 0.2)"
                    height={15}
                    width={40}
                    className='react-switch darkmode-switch'
                />
            </JSONViewer>
        </SplitPane>
    );
};

export default CalculatorWrapper;
