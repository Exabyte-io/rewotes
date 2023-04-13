import React, { useState, useCallback } from 'react';
import { applyEdgeChanges, applyNodeChanges } from 'reactflow';
import createNode from '../utils/createNode';
import useUpdateOutputNodes from './useUpdateOutputNodes';

const useLocalFlowData = (initial = {}) => {
    const [reactFlowInstance, setReactFlowInstance] = useState(null);
    const [nodes, setNodes] = useState(initial.nodes || []);
    const [edges, setEdges] = useState(initial.edges || []);

    const updateOutputNodes = useUpdateOutputNodes(nodes, edges, setNodes);

    // handle nodes and edges changes
    const onNodesChange = useCallback(
        (changes) => setNodes((els) => applyNodeChanges(changes, els)),
        []
    );

    const onEdgesChange = useCallback(
        (changes) => setEdges((els) => applyEdgeChanges(changes, els)),
        []
    );

    // Function to add a new node
    const addNode = (type, position) => {
        const newNode = createNode(type, position, undefined, handleNodeChange);
        setNodes((nodes) => [...nodes, newNode]);
    };

    const updateReactFlowInstance = (instance) => {
        setReactFlowInstance(instance);
    };

    const clearFlowchart = () => {
        setNodes([]);
        setEdges([]);
    };

    const loadFlowchart = (flow) => {
        const newNodes = attachOnConnect(flow.nodes);
        setNodes(newNodes);
        setEdges(flow.edges);
    };

    const handleNodeChange = (id, value) => {
        setNodes((ns) => {
            const newNodes = ns.map((n) =>
                n.id === id ? { ...n, data: { ...n.data, value } } : n
            );
            return newNodes;
        });
    };

    const attachOnConnect = (nodes) => {
        return nodes.map((node) => ({
            ...node,
            data: {
                ...node.data,
                onChange: (newValue) => {
                    handleNodeChange(node.id, newValue);
                },
            },
        }));
    };

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

    return {
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
    };
};

export default useLocalFlowData;
