import React, { useEffect, useMemo } from 'react';
import ReactFlow, { Background, Controls, MiniMap } from 'reactflow';

import nodeTypesConfig from '../customNodes/nodeTypes';
import usePrevious from '../../hooks/usePrevious';
import getNodeColor from '../../utils/getNodeColor';
import { useDarkMode } from '../reusable/DarkModeContext';

const FlowchartCanvas = ({
    setReactFlowInstance,
    nodes,
    edges,
    onNodesChange,
    onEdgesChange,
    handleConnect,
    handleDragOver,
    handleDrop,
    updateOutputNodes,
    // isDarkMode,
}) => {
    const { isDarkMode } = useDarkMode();
    // store previous nodes and edges as state
    const prevNodes = usePrevious(nodes);
    const prevEdges = usePrevious(edges);

    // Define custom node types
    const nodeTypes = useMemo(() => {
        return nodeTypesConfig;
    }, []);

    useEffect(() => {
        window.handleConnect = handleConnect;
        window.edges = edges;
        return () => {
            window.handleConnect = null;
            window.edges = null;
        };
    }, [handleConnect, edges]);

    useEffect(() => {
        // TODO: use lodash to evaluate equality?
        if (
            JSON.stringify(prevNodes) !== JSON.stringify(nodes) ||
            JSON.stringify(prevEdges) !== JSON.stringify(edges)
        ) {
            updateOutputNodes();
        }
    }, [nodes, edges, prevNodes, prevEdges, updateOutputNodes]);

    return (
        <ReactFlow
            nodes={nodes}
            edges={edges}
            onNodesChange={onNodesChange}
            onEdgesChange={onEdgesChange}
            onConnect={handleConnect}
            onInit={setReactFlowInstance}
            nodeTypes={nodeTypes}
            snapToGrid
            onDrop={handleDrop}
            onDragOver={handleDragOver}
            proOptions={{ hideAttribution: true }}
            className={isDarkMode ? "dark-mode" : ""}
        >
            <MiniMap nodeColor={getNodeColor} />
            <Controls />
            <Background className={isDarkMode ? "dark-mode" : ""}
                gap={16}
                color={
                    isDarkMode
                        ? 'rgba(255, 255, 255, 1)'
                        : 'rgba(100, 100, 100, 1)'
                }
            />
        </ReactFlow>
    );
};

export default FlowchartCanvas;
