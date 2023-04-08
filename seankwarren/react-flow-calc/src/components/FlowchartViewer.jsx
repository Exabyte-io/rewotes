import React, { useCallback, useMemo } from 'react'
import ReactFlow, {
    Controls,
    Background,
    MiniMap,
    applyNodeChanges,
    applyEdgeChanges,
} from 'reactflow'
import InputNode from './InputNode';
import OperationNode from './OperationNode';
import ComparisonNode from './ComparisonNode';
import OutputNode from './OutputNode';

const FlowchartViewer = ({ 
    nodes, 
    edges, 
    setNodes,
    setEdges, 
    addNode, 
    draggedNodeType,
    setDraggedNodeType, 
}) => {

    const onNodesChange = useCallback(
        (changes) => setNodes((els) => applyNodeChanges(changes, els)),
        []
    );

    const onEdgesChange = useCallback(
        (changes) => setEdges((els) => applyEdgeChanges(changes, els)),
        []
    );

    const handleDrop = (e) => {
        console.log('dropped', draggedNodeType);
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

    const nodeTypes = useMemo(() => {
        return {
            inputNode: (props) => <InputNode {...props}/>,
            operationNode: (props) => <OperationNode {...props}/>,
            outputNode: (props) => <OutputNode {...props} />,
            comparisonNode: (props) => <ComparisonNode {...props}/>
        };
    }, []);

    return (
        <ReactFlow
            nodes={nodes}
            edges={edges}
            onNodesChange={onNodesChange}
            onEdgesChange={onEdgesChange}
            nodeTypes={nodeTypes}
            snapToGrid
            onDrop={handleDrop}
            onDragOver={handleDragOver}
        >
            <MiniMap />
            <Controls />
            <Background />
        </ReactFlow>
    )
}

export default FlowchartViewer