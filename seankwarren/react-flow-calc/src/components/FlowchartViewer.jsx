import React, { useCallback, useEffect, useMemo } from 'react';
import ReactFlow, {
    Controls,
    Background,
    MiniMap,
    applyNodeChanges,
    applyEdgeChanges,
} from 'reactflow';
import InputNode from './InputNode';
import OperationNode from './OperationNode';
import OutputNode from './OutputNode';
import ComparisonNode from './ComparisonNode';
import calculate from '../utils/calculate';
import usePrevious from '../hooks/usePrevious';

const FlowchartViewer = ({
    nodes,
    edges,
    setNodes,
    setEdges,
    addNode,
    draggedNodeType,
    setDraggedNodeType,
}) => {
    const prevNodes = usePrevious(nodes);
    const prevEdges = usePrevious(edges);

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
        // If a duplicate edge is found, do not add the new edge
        if (existingEdge) {
            console.log('Duplicate edge detected, not adding the new edge');
            return;
        }

        // Create a new edge with an ID, source node ID, target node ID, and arrowhead type
        const newEdge = {
            id: `${source}-${target}`,
            source: source,
            sourceHandle: sourceHandle,
            target: target,
            targetHandle: targetHandle,
            arrowHeadType: 'arrowclosed',
        };

        const newEdges = [...edges, newEdge];
        setEdges(newEdges);

        // If the target node is an output node, calculate its value and update its state
        updateOutputNodes();
    };

    const updateOutputNodes = (currentEdges) => {
        setNodes((currentNodes) => {
            const outputNodes = currentNodes.filter((node) => {
                return node.type === 'outputNode';
            });
            const newNodes = currentNodes.map((node) => {

                if (node.type !== 'outputNode') return node;
                
                const connectedEdge = edges.find(
                    (edge) => edge.target === node.id
                );
                if (connectedEdge) {
                    const newValue = calculate(
                        nodes,
                        edges,
                        connectedEdge.sourceHandle
                    );
                    return { ...node, data: { ...node.data, value: newValue } };
                } else {
                    return node;
                }
            });
            return newNodes;
        });
    }

    // Define custom node types
    const nodeTypes = useMemo(() => {
        return {
            inputNode: (props) => (
                <InputNode 
                    {...props} 
                />
            ),
            operationNode: (props) => (
                <OperationNode
                    {...props}
                />
            ),
            comparisonNode: (props) => (
                <ComparisonNode
                    {...props}
                />
            ),
            outputNode: (props) => <OutputNode {...props} />,
        };
    }, []);

    useEffect(() => {
        // TODO: use lodash to evaluate equality?
        if (
            JSON.stringify(prevNodes) !== JSON.stringify(nodes) ||
            JSON.stringify(prevEdges) !== JSON.stringify(edges)
        ) {
            updateOutputNodes();
        }
    }, [nodes, edges, prevNodes, prevEdges]);

    const handleDrop = (e) => {
        e.preventDefault();
        e.stopPropagation();

        if (!draggedNodeType) return;

        const reactFlowBounds = e.currentTarget.getBoundingClientRect();
        const position = {
            x: e.clientX - reactFlowBounds.left,
            y: e.clientY - reactFlowBounds.top,
        };

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
            nodeTypes={nodeTypes}
            snapToGrid
            onDrop={handleDrop}
            onDragOver={handleDragOver}
        >
            <MiniMap />
            <Controls />
            <Background />
        </ReactFlow>
    );
};

export default FlowchartViewer;
