import {Node, Edge} from 'react-flow-renderer'

export const SET_NODES_DATA = 'SET_NODES_DATA'
export const SET_EDGES_DATA = 'SET_EDGES_DATA'

export const setNodesData = (payload: Node[]): Action => ({
    type: SET_NODES_DATA,
    payload: payload
})

export const setEdgesData = (payload: Edge[]): Action => ({
    type: SET_EDGES_DATA,
    payload: payload
})

export interface Action {
    type: string,
    payload: any
}