import React, { useCallback, useEffect, useMemo, useState } from 'react';
import ReactFlow, {
    Background,
    Controls,
    MiniMap,
    applyEdgeChanges,
    applyNodeChanges,
} from 'reactflow';
import usePrevious from '../../hooks/usePrevious';
import useUpdateOutputNodes from '../../hooks/useUpdateOutputNodes';
import getNodeColor from '../../utils/getNodeColor';
import nodeTypesConfig from '../customNodes/nodeTypes';

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
    isDarkMode,
}) => {
    // store previous nodes and edges as state
    const prevNodes = usePrevious(nodes);
    const prevEdges = usePrevious(edges);
    
    // const updateOutputNodes = useUpdateOutputNodes(nodes, edges, setNodes);

    // Darkmode style toggling
    const reactFlowStyle = {
        backgroundColor: isDarkMode
            ? 'rgba(40, 40, 40, 1)'
            : 'rgba(255, 255, 255, 1)',
    };

    // Define custom node types
    const nodeTypes = useMemo(() => {
        return nodeTypesConfig
    }, []);

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
            style={reactFlowStyle}
            proOptions={{ hideAttribution: true }}
        >
            <MiniMap nodeColor={getNodeColor} />
            <Controls />
            <Background
                gap={16}
                color={
                    isDarkMode
                        ? 'rgba(255, 255, 255, 1)'
                        : 'rgba(100, 100, 100, 1)'
                }
                isDarkMode={isDarkMode}
            />
        </ReactFlow>
    );
};

export default FlowchartCanvas;
