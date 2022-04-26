import {Node, Edge} from 'react-flow-renderer'

export const SET_NODES_DATA = 'SET_NODES_DATA'
export const SET_EDGES_DATA = 'SET_EDGES_DATA'
export const UPDATE_DATA_BY_ID = 'UPDATE_NODE_BY_ID'
export const UPDATE_NODES = 'UPDATE_NODES'

export const setNodesData = (payload: Node[]): Action => ({
    type: SET_NODES_DATA,
    payload: payload
})

export const updateDataById = (payload: {id: string, data: string}): Action => ({
    type: UPDATE_DATA_BY_ID,
    payload: payload
})

export const updateNodes = (payload: Node[]): Action => ({
    type: UPDATE_NODES,
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