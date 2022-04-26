import {Node, Edge} from 'react-flow-renderer'
import {Action, SET_EDGES_DATA, SET_NODES_DATA, UPDATE_DATA_BY_ID, UPDATE_NODES} from "./actions";

export const nodesReducer = (nodes: Node[] = [], action: Action): Node[] => {
    switch (action.type) {
        case SET_NODES_DATA:
            return action.payload
        case UPDATE_DATA_BY_ID:
            const index = nodes.findIndex(item => item.id === action.payload.id)
            if (index >= 0) {
                nodes[index].data = action.payload.data
                return [...nodes]
            }
            return nodes
        case UPDATE_NODES:
            return nodes.map((node: Node) => {
                if (action.payload.map((item: Node) => item.id).includes(node.id)) {
                    return action.payload.find((item: Node) => item.id === node.id)
                }
                return node
            })
    }
    return nodes
}

export const edgesReducer = (edges: Edge[] = [], action: Action): Edge[] => action.type === SET_EDGES_DATA ? action.payload : edges