import {Node, Edge} from 'react-flow-renderer'
import {Action, SET_EDGES_DATA, SET_NODES_DATA} from "./actions";

export const nodesReducer = (nodes: Node[] = [], action: Action): Node[] => action.type === SET_NODES_DATA ? action.payload : nodes

export const edgesReducer = (edges: Edge[] = [], action: Action): Edge[] => action.type === SET_EDGES_DATA ? action.payload : edges