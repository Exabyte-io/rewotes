import React, {useState, useCallback} from "react";
import ReactFlow, {Node, Edge, NodeChange, EdgeChange, applyNodeChanges, applyEdgeChanges, addEdge, Connection, MiniMap, Background, Handle, Position, NodeTypes } from 'react-flow-renderer';
import nodeTypes from "../components/customNodes/CustomNodeTypes";

type Props = {
    initialNodes: Node[],
    initialEdges: Edge[]
}

const FlowChart: React.FC<Props> = (props: Props) => {
    const {initialNodes, initialEdges} = props
    const [nodes, setNodes] = useState<Node[]>(initialNodes);
    const [edges, setEdges] = useState<Edge[]>(initialEdges);

    const onNodesChange = useCallback(
        (changes: NodeChange[]) => setNodes((nds) => applyNodeChanges(changes, nds)),
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

    return <ReactFlow
        nodeTypes={nodeTypes}
        nodes={nodes}
        edges={edges}
        onNodesChange={onNodesChange}
        onEdgesChange={onEdgesChange}
        onConnect={onConnect}
        fitView >
        <Background />
    </ReactFlow>
}

export default FlowChart