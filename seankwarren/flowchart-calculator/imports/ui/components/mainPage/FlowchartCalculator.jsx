import React, { useState } from 'react';
import Switch from 'react-switch';
import 'reactflow/dist/style.css';
import SplitPane from 'split-pane-react';
import 'split-pane-react/esm/themes/default.css';
import createNode from '../../utils/createNode';
import startingNode from '../../utils/startingNode';
import FlowchartCanvas from './FlowchartCanvas';
import JSONViewer from './JSONViewer';
import NodeButtons from './NodeButtons';

const FlowchartCalculator = () => {
    // Set up state for nodes and edges
    const [nodes, setNodes] = useState([startingNode]);
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

    const handleNodeChange = (id, value) => {
        setNodes((ns) => {
            const newNodes = ns.map((n) => 
                n.id === id
                    ? { ...n, data: { ...n.data, value } }
                    : n
            );
            return newNodes;
        });
    };

    // Function to add a new node
    const addNode = (type, position) => {
        const newNode = createNode(type, position, undefined, handleNodeChange);
        setNodes((nodes) => [...nodes, newNode]);
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
                <FlowchartCanvas
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

export default FlowchartCalculator;
