import React, {useEffect} from "react";
import ReactFlow, {
    addEdge,
    applyEdgeChanges,
    applyNodeChanges,
    Background,
    Connection,
    Edge,
    EdgeChange,
    Node,
    NodeChange
} from 'react-flow-renderer';
import {useDispatch, useSelector} from "react-redux";
import nodeTypes, {operatorNodes} from "../components/customNodes/CustomNodeTypes";
import {ReducersType} from "../redux/store";
import {setEdgesData, setNodesData} from "../redux/actions";


const FlowChart: React.FC = () => {
    const dispatch = useDispatch()
    const nodes: Node[] = useSelector((state: ReducersType) => state.nodes)
    const edges: Edge[] = useSelector((state: ReducersType) => state.edges)

    useEffect(() => {
        localStorage.setItem('nodes', JSON.stringify(nodes))
        localStorage.setItem('edges', JSON.stringify(edges))
    }, [nodes, edges])

    const onNodesChange = (changes: NodeChange[]) => {
        dispatch(setNodesData(applyNodeChanges(changes, nodes)))
    }

    const onEdgesChange = (changes: EdgeChange[]) => {
        dispatch(setEdgesData(applyEdgeChanges(changes, edges)))
    }
    const onConnect = (connection: Connection) => {
        const source: Node | undefined = nodes.find((item: Node) => item.id === connection.source)
        if (source && source.type && operatorNodes.includes(source.type)) {
            const edgesOfSource: Edge[] = edges.filter((item: Edge) => item.source === source.id)
            if (!edgesOfSource.find(item => item.sourceHandle === connection.sourceHandle)) {
                dispatch(setEdgesData(addEdge(connection, edges)))
            }
        } else {
            dispatch(setEdgesData(addEdge(connection, edges)))
        }
    }

    return <ReactFlow
        nodeTypes={nodeTypes}
        nodes={nodes}
        edges={edges}
        onNodesChange={onNodesChange}
        onEdgesChange={onEdgesChange}
        onConnect={onConnect}
        fitView>
        <Background/>
    </ReactFlow>
}

export default FlowChart