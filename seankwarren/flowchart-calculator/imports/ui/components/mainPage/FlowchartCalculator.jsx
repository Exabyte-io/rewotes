import React, { useEffect } from 'react';

import 'reactflow/dist/style.css';
import 'split-pane-react/esm/themes/default.css';

import FlowchartCanvas from './FlowchartCanvas';
import JSONViewer from './JSONViewer';
import NodeButtons from './NodeButtons';
import useLocalFlowData from '../../hooks/useLocalFlowData';
import useRemoteFlowData from '../../hooks/useRemoteFlowData';
import useDraggable from '../../hooks/useDraggable';
import ResizablePane from '../reusable/ResizablePane';
import startingNode from '../../utils/startingNode';
import { DarkModeProvider } from '../reusable/DarkModeContext';
import DarkModeSwitch from '../reusable/DarkModeSwitch';

const FlowchartCalculator = () => {
    // state for nodes and edges
    const {
        reactFlowInstance,
        updateReactFlowInstance,
        clearFlowchart,
        loadFlowchart,
        nodes,
        onNodesChange,
        addNode,
        edges,
        onEdgesChange,
        handleConnect,
        updateOutputNodes,
    } = useLocalFlowData({ nodes: [startingNode] });

    // state for flows on db
    const {
        fetchedFlows,
        saveFlow,
        fetchFlows,
        clearFlows,
        flowName,
        updateFlowName,
    } = useRemoteFlowData(nodes, edges);

    // Drag and drop handling
    const { handleDragStart, handleDragOver, handleDrop } = useDraggable(
        reactFlowInstance,
        addNode
    );

    useEffect(() => {
        fetchFlows();
    }, []);

    return (
        <DarkModeProvider>
            <ResizablePane>
                <div
                    className='flowchart-container'
                    style={{ height: '100vh', width: '100%' }}
                >
                    <NodeButtons
                        handleDragStart={handleDragStart}
                        clearFlows={clearFlows}
                        clearFlowchart={clearFlowchart}
                    />
                    <FlowchartCanvas
                        setReactFlowInstance={updateReactFlowInstance}
                        nodes={nodes}
                        edges={edges}
                        onNodesChange={onNodesChange}
                        onEdgesChange={onEdgesChange}
                        addNode={addNode}
                        handleConnect={handleConnect}
                        handleDragOver={handleDragOver}
                        handleDrop={handleDrop}
                        updateOutputNodes={updateOutputNodes}
                    />
                </div>
                <JSONViewer
                    nodes={nodes}
                    edges={edges}
                    flows={fetchedFlows}
                    loadFlowchart={loadFlowchart}
                    onSave={saveFlow}
                    flowName={flowName}
                    updateFlowName={updateFlowName}
                >
                    <DarkModeSwitch />
                </JSONViewer>
            </ResizablePane>
        </DarkModeProvider>
    );
};

export default FlowchartCalculator;
