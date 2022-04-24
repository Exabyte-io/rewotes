import React, {useEffect, useState} from "react";
import {Edge, Node} from "react-flow-renderer";
import {store} from "../redux/store";
import {useDispatch} from "react-redux";
import {updateDataById} from "../redux/actions";

const RPN: React.FC = () => {
    const dispatch = useDispatch()

    const [nodes, setNodes] = useState<Node[]>([])
    const [edges, setEdges] = useState<Edge[]>([])

    const [operators, setOperators] = useState<Node[]>([])
    const [IONodes, setIONodes] = useState<Node[]>([])
    const [connectionsOfOperators, setConnectionsOfOperators] = useState<Edge[]>([])

    useEffect(() => {
        const interval = setInterval(() => {
            setNodes(store.getState().nodes)
            setEdges(store.getState().edges)
        }, 1000)

        return () => clearInterval(interval)
    }, [])

    useEffect(() => {
        setOperators(getAllOperators(nodes))
        setIONodes(getIONodes(nodes))
    }, [nodes])

    useEffect(() => {
        setConnectionsOfOperators(getConnectionsOfOperators(operators, edges))
    }, [operators, edges])

    useEffect(() => {
        getCalculatedArrayOfNodes(operators, IONodes, connectionsOfOperators).forEach((item: Node) => {
            dispatch(updateDataById({id: item.id, data: item.data}))
        })
    }, [operators, IONodes, connectionsOfOperators])

    function getAllOperators(nodes: Node[]): Node[] {
        return nodes.filter(item => ['plus', 'minus', 'divide', 'multiply'].includes(item.type as string))
    }

    function getIONodes(nodes: Node[]): Node[] {
        return nodes.filter(item => item.type === 'io')
    }

    function getConnectionsOfOperators(nodes: Node[], edges: Edge[]): Edge[] {
        const ids: string[] = nodes.reduce((res: string[], item) => {
            res.push(item.id)
            return res
        }, [])
        return edges.filter(item => ids.includes(item.source))
    }

    function getCalculatedArrayOfNodes(operators: Node[], IONodes: Node[], connectionsOfOperators: Edge[]): Node[] {
        const arrayOfRelations: Node[][] = operators.reduce((res: Node[][], node) => {
            const targetNodes: string[] = connectionsOfOperators.filter(item => item.source === node.id).map(item => item.target)
            const targetIOs = IONodes.filter(item => targetNodes.includes(item.id))
            const targetOperators = operators.filter(item => targetNodes.includes(item.id))
            const sortedTargetsIds: string[] = connectionsOfOperators.filter(item => item.source === node.id).sort((a, b) => `${a.sourceHandle}` > `${b.sourceHandle}` ? 1 : -1).map(item => item.target)
            const children: Node[] = sortedTargetsIds.reduce((res: Node[], id) => {
                const node: Node = [...operators, ...IONodes].find(item => item.id === id) as Node
                if (node) res.push(node)
                return res
            }, [])
            res.push([node, ...children])
            return res
        }, [])

        return [...calculateNodes(arrayOfRelations)]
    }

    function calculateNodes(arr: Node[][]): Node[] {
        let recall: boolean = false
        arr.forEach((nodes: Node[]) => {
            let isCalculated: boolean = true
            nodes.forEach((node: Node, index: number) => {
                if (index !== 0 && ['plus', 'minus', 'divide', 'multiply'].includes(node.type as string) && typeof node.data !== "number") {
                    const isChildOperatorHasAChildren = arr.reduce((res: boolean, item) => {
                        if (item[0].id === node.id && item.length > 1) res = true
                        return res
                    }, false)
                    if (isChildOperatorHasAChildren) {
                        isCalculated = false
                        recall = true
                    }
                }
            })
            if (isCalculated) {
                nodes[0].data = calculateValues(nodes)
            }
        })
        if (recall) {
            calculateNodes(arr)
        } else {
            return arr.reduce((res: Node[], item) => {
                res.push(...item)
                return res
            }, [])
        }
        return []
    }

    function calculateValues(nodes: Node[]): number {
        const typeOfOperator = nodes[0].type
        return  nodes.reduce((res: number, item: Node, index: number) => {
            const data = item.data === '' ? 0 : + item.data
            if (index === 1) res = data
            if (index !== 0 && index !== 1) {
                switch (typeOfOperator) {
                    case 'plus': res += data
                        return res
                    case 'minus': res = res - data
                        return res
                    case 'divide': res = res / data
                        return res
                    case 'multiply': res *= data
                        return res
                }
            }
            return res
        }, 0)
    }



    return <></>
}

export default RPN