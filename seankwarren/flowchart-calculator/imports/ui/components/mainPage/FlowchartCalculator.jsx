import React, { useState, useEffect } from 'react';
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
    // state for nodes and edges
    const [nodes, setNodes] = useState([startingNode]);
    const [edges, setEdges] = useState([]);
    const [reactFlowInstance, setReactFlowInstance] = useState(null);

    // state for fetching saved flows from db
    const [fetchedFlows, setFetchedFlows] = useState([]);

    // state for user input flow name
    const [flowName, setFlowName] = useState('');


    const saveFlow = () => {
        // Check if a flow with the same name already exists
        const flowWithNameExists = fetchedFlows.some(flow => flow.name === flowName);
        
        if (flowWithNameExists) {
            alert('A flow with this name already exists. Please choose a different name.');
        } else {
            Meteor.call('saveFlow', {nodes, edges, name: flowName}, (error) => {
                if (error) {
                    console.log(error.reason);
                } else {
                    console.log('Data saved to MongoDB');
                    setFlowName(''); // Clear the flow name input field after a successful save
                    fetchFlows(); // Fetch the flows again to update the dropdown
                }
            });
        }
    };

    const fetchFlows = () => {
        Meteor.call('fetchFlows', (error, result) => {
            if (error) {
                console.error("Error fetching flows:", error);
            } else {
                setFetchedFlows(result);
            }
        });
    };

    const clearFlows = () => {
        Meteor.call('clearFlows', (error) => {
            if (error) {
                console.error("Error clearing flows:", error);
            } else {
                console.log("Flows collection cleared");
            }
        });
    };


    const clearFlowchart = () => {
        console.log("clearing flowchart");
        setReactFlowInstance(null);
        setNodes([]);
        setEdges([]);
    }

    const loadFlow = (flow) => {
        const newNodes = attachOnConnect(flow.nodes);
        setNodes(newNodes);
        setEdges(flow.edges);
    };
      

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

    const attachOnConnect = (nodes) => {
        return nodes.map((node) => ({
            ...node,
            data: { ...node.data, onChange: (newValue) => { 
                handleNodeChange(node.id, newValue) 
            } },
        }));
    };

    useEffect(() => {
        fetchFlows();
    }, []);

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
                    clearFlows={clearFlows}
                    clearFlowchart={clearFlowchart}
                />
                <FlowchartCanvas
                    nodes={nodes}
                    edges={edges}
                    setNodes={setNodes}
                    setEdges={setEdges}
                    addNode={addNode}
                    reactFlowInstance={reactFlowInstance}
                    setReactFlowInstance={setReactFlowInstance}
                    draggedNodeType={draggedNodeType}
                    setDraggedNodeType={setDraggedNodeType}
                    isDarkMode={isDarkMode}
                />
            </div>
            <JSONViewer 
                nodes={nodes} 
                edges={edges} 
                isDarkMode={isDarkMode}
                flows={fetchedFlows}
                loadFlow={loadFlow}
                saveFlow={saveFlow}
                onSave={saveFlow}
                flowName={flowName}
                setFlowName={setFlowName}
            >
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
