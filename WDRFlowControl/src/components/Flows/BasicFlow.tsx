import ReactFlow, {
    addEdge,
    applyEdgeChanges,
    applyNodeChanges,
    Background,
    Connection, Controls,
    Edge,
    EdgeChange, FitViewOptions,
    Node,
    NodeChange
} from "reactflow";
import React, { useCallback, useMemo, useState } from "react";
import NumberInput from "../customNodes/NumberInput";
import OperationInput from "../customNodes/OperationInput";
import TestNode from "../customNodes/TestNode";
import AdjustmentInput from "../customNodes/AdjustmentInput";

export const BasicFlow = ({nodes, edges, setNodes, setEdges}: {
    nodes:Node[]
    setNodes: (p: (nds: Node[]) => Node[]) => void
    setEdges: (p: (nds: Edge[]) => Edge[]) => void
    edges: Edge[]
}): JSX.Element => {
    const onNodesChange = useCallback(
        (changes: NodeChange[]) => {
            setNodes((nds) => applyNodeChanges(changes, nds))
        },
        [setNodes]
    );
    const onEdgesChange = useCallback(
        (changes: EdgeChange[]) => setEdges((eds) => applyEdgeChanges(changes, eds)),
        [setEdges]
    );
    const onConnect = useCallback(
        (connection: Connection) => setEdges((eds) => addEdge(connection, eds)),
        [setEdges]
    );

    const nodeTypes = useMemo(() => ({
        numberInput: NumberInput,
        operationInput: OperationInput,
        testNode: TestNode,
        adjustmentInput: AdjustmentInput,
    }), [])


    const fitViewOptions: FitViewOptions = {
        padding: 0.2,
    };

    return (
        <ReactFlow
            nodes={nodes}
            edges={edges}
            onNodesChange={onNodesChange}
            onEdgesChange={onEdgesChange}
            onConnect={onConnect}
            fitView
            fitViewOptions={fitViewOptions}
            nodeTypes={nodeTypes}
        >
            <Background/>
            <Controls/>
        </ReactFlow>
    )
}
