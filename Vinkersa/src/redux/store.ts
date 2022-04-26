import { createStore, combineReducers} from 'redux'
import {Node, Edge} from 'react-flow-renderer'
import {edgesReducer, nodesReducer} from "./reducers";

export interface ReducersType {
    nodes: Node[],
    edges: Edge[]
}

const allReducers = combineReducers<ReducersType>({
    nodes: nodesReducer,
    edges: edgesReducer
})

export const store = createStore(allReducers)
