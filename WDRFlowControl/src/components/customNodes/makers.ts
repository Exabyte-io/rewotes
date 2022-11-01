// as functions so we could create "Add Input" buttons
import { Node } from "reactflow";

export const makeNumberInput = ({ id, data, position }: Node): Node => {
    return ({ type: 'numberInput', id, data, position })
}

export const makeOperationInput = ({ id, data, position }: Node): Node => {
    return ({ type: 'operationInput', id, data, position })
}

export const makeTestNode = (props: Node): Node => {
    const { id, data, position } = props
    return ({ type: 'testNode', id, data, position })
}

export const makeAdjustmentInput = (props: Node): Node => {
    const { id, data, position } = props
    return ({ type: 'adjustmentInput', id, data, position })
}