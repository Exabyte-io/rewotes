import React, {useState, useCallback} from "react";
import ReactFlow, {Node, Edge, NodeChange, EdgeChange, applyNodeChanges, applyEdgeChanges, addEdge, Connection, MiniMap, Background, Handle, Position, NodeTypes } from 'react-flow-renderer';
import nodeTypes from "../components/customNodes/CustomNodeTypes";

const initialNodes: Node[] = [
    {
        id: '1',
        type: 'terminal',
        data: { label: 'Input Node1' },
        position: { x: 250, y: 25 },
    },

    {
        id: '2',
        data: { label: <div>Default Node</div> },
        position: { x: 100, y: 125 },
    },
    {
        id: '3',
        type: 'output',
        data: { label: 'Output Node' },
        position: { x: 250, y: 250 },
    }
];

const initialEdges: Edge[] = [
    { id: 'e1-2', source: '1', target: '2' },
    { id: 'e2-3', source: '2', target: '3', animated: true },
];


const FlowChart: React.FC = () => {
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
        <MiniMap />
    </ReactFlow>
}

export default FlowChart