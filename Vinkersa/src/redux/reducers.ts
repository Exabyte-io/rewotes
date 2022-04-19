import {Node, Edge} from 'react-flow-renderer'
import {Action, SET_EDGES_DATA, SET_NODES_DATA, UPDATE_DATA_BY_ID} from "./actions";

export const nodesReducer = (nodes: Node[] = [], action: Action): Node[] => {
    switch (action.type) {
        case SET_NODES_DATA: return action.payload
        case UPDATE_DATA_BY_ID:
            const index = nodes.findIndex(item => item.id === action.payload.id)
            if (index >= 0) {
                nodes[index].data = action.payload.data
                return [...nodes]
            }
            return nodes
    }
    return  nodes
}

export const edgesReducer = (edges: Edge[] = [], action: Action): Edge[] => action.type === SET_EDGES_DATA ? action.payload : edges