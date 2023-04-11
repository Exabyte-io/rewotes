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
    nodes,
    edges,
    setNodes,
    setEdges,
    addNode,
    draggedNodeType,
    setDraggedNodeType,
    isDarkMode,
}) => {
    const [reactFlowInstance, setReactFlowInstance] = useState(null);

    // store previous nodes and edges as state
    const prevNodes = usePrevious(nodes);
    const prevEdges = usePrevious(edges);
    
    const updateOutputNodes = useUpdateOutputNodes(nodes, edges, setNodes);

    // Darkmode style toggling
    const reactFlowStyle = {
        backgroundColor: isDarkMode
            ? 'rgba(40, 40, 40, 1)'
            : 'rgba(255, 255, 255, 1)',
    };

    // handle nodes and edges changes
    const onNodesChange = useCallback(
        (changes) => setNodes((els) => applyNodeChanges(changes, els)),
        []
    );

    const onEdgesChange = useCallback(
        (changes) => setEdges((els) => applyEdgeChanges(changes, els)),
        []
    );

    // Function to handle connecting nodes
    const handleConnect = (params) => {
        const { source, sourceHandle, target, targetHandle } = params;

        // Check if an edge with the same source and target already exists
        const existingEdge = edges.find(
            (edge) =>
                edge.source === sourceHandle && edge.target === targetHandle
        );
        if (existingEdge) return;

        // Check if there is already an edge connected to the target handle
        const existingEdgeWithSameTargetHandle = edges.find(
            (edge) => edge.targetHandle === targetHandle
        );
        if (existingEdgeWithSameTargetHandle) return;

        // Create a new edge with an ID, source node ID, target node ID, and arrowhead type
        const newEdge = {
            id: `${source}-${target}`,
            source: source,
            sourceHandle: sourceHandle,
            target: target,
            targetHandle: targetHandle,
            arrowHeadType: 'arrowclosed',
            animated: true,
            interactionWidth: 40,
        };

        const newEdges = [...edges, newEdge];
        setEdges(newEdges);

        // If the target node is an output node, calculate its value and update its state
        updateOutputNodes();
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

    const handleDrop = (e) => {
        e.preventDefault();
        e.stopPropagation();

        if (!draggedNodeType) return;

        const reactFlowBounds = e.currentTarget.getBoundingClientRect();
        const position = reactFlowInstance.project({
            x: e.clientX - reactFlowBounds.left,
            y: e.clientY - reactFlowBounds.top,
        });

        addNode(draggedNodeType, position);
        setDraggedNodeType(null);
    };

    const handleDragOver = (e) => {
        e.preventDefault();
        e.stopPropagation();
    };

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
