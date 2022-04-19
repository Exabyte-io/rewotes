import React, {useEffect} from "react";
import ReactFlow, {
    addEdge,
    applyEdgeChanges,
    applyNodeChanges,
    Background,
    Connection,
    ConnectionMode,
    Edge,
    EdgeChange,
    Node,
    NodeChange
} from 'react-flow-renderer';
import {useDispatch, useSelector} from "react-redux";
import nodeTypes from "../components/customNodes/CustomNodeTypes";
import {ReducersType} from "../redux/store";
import {setEdgesData, setNodesData} from "../redux/actions";


const FlowChart: React.FC = () => {
    const dispatch = useDispatch()
    const nodes: Node[] = useSelector((state: ReducersType) => state.nodes)
    const edges: Edge[] = useSelector((state: ReducersType) => state.edges)

    useEffect(() => {
        if (nodes.length > 0) localStorage.setItem('nodes', JSON.stringify(nodes))
        if (edges.length > 0) localStorage.setItem('edges', JSON.stringify(edges))
    }, [nodes, edges])

    const onNodesChange = (changes: NodeChange[]) => {
        dispatch(setNodesData(applyNodeChanges(changes, nodes)))
    }

    const onEdgesChange = (changes: EdgeChange[]) => {
        dispatch(setEdgesData(applyEdgeChanges(changes, edges)))
    }
    const onConnect = (connection: Connection) => {
        dispatch(setEdgesData(addEdge(connection, edges)))
    }

    return <ReactFlow
        nodeTypes={nodeTypes}
        nodes={nodes}
        edges={edges}
        onNodesChange={onNodesChange}
        onEdgesChange={onEdgesChange}
        onConnect={onConnect}
        connectionMode={ConnectionMode.Loose}
        fitView>
        <Background/>
    </ReactFlow>
}

export default FlowChart