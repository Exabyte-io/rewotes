import React, { useCallback, useEffect, useMemo } from 'react';
import ReactFlow, {
    Controls,
    Background,
    MiniMap,
    applyNodeChanges,
    applyEdgeChanges,
} from 'reactflow';
import InputNode from './customNodes/InputNode';
import BinaryNode from './customNodes/BinaryNode';
import UnaryNode from './customNodes/UnaryNode';
import OutputNode from './customNodes/OutputNode';
import ComparisonNode from './customNodes/ComparisonNode';
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
    isDarkMode
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
            binaryNode: (props) => (
                <BinaryNode
                    {...props}
                />
            ),
            unaryNode: (props) => (
                <UnaryNode
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

    const reactFlowStyle = {
        backgroundColor: isDarkMode ? 'rgba(40, 40, 40, 1)' : 'rgba(255, 255, 255, 1)',
    };

    const customNodeColor = (node) => {
        switch (node.type) {
            case 'inputNode':
                return 'var(--input-color)';
            case 'binaryNode':
                return 'var(--binary-color)';
            case 'unaryNode':
                return 'var(--unary-color)';
            case 'comparisonNode':
                return 'var(--comparison-color)';
            case 'outputNode':
                return 'var(--output-color)';
            default:
                return '#eee';
        }
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
            style={reactFlowStyle}
        >
            <MiniMap nodeColor={customNodeColor} />
            <Controls />
            <Background 
                gap={16}
                color={isDarkMode ? 'rgba(255, 255, 255, 1)' : 'rgba(100, 100, 100, 1)'} 
                isDarkMode={isDarkMode}
            />
        </ReactFlow>
    );
};

export default FlowchartViewer;
