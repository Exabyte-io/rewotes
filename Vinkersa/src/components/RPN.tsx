import React, {useEffect, useState} from "react";
import {Edge, Node} from "react-flow-renderer";
import {ReducersType, store} from "../redux/store";
import {useDispatch, useSelector} from "react-redux";
import {updateDataById} from "../redux/actions";
import {operatorNodes} from "./customNodes/CustomNodeTypes";

const RPN: React.FC = () => {
    const dispatch = useDispatch()

    const nodes: Node[] = useSelector((state: ReducersType) => state.nodes)
    const edges: Edge[] = useSelector((state: ReducersType) => state.edges)

    const [operators, setOperators] = useState<Node[]>([])
    const [IONodes, setIONodes] = useState<Node[]>([])

    const [prev, setPrev] = useState<{state: string, isChanged: boolean}>({state: '', isChanged: false})

    useEffect(() => {
        setOperators(getAllOperators(nodes))
        setIONodes(getIONodes(nodes, edges))
    }, [nodes, edges])

    useEffect(() => {
        const data: string = JSON.stringify([operators.map(item => item.data), IONodes.map(item => item.data)])
        if (data !== prev.state) {
            setPrev({state: data, isChanged: true})
        }
    }, [operators, IONodes])

    useEffect(() => {
        if (prev.isChanged) {
            getCalculatedArrayOfNodes(operators, IONodes, edges).forEach((item: Node) => {
                dispatch(updateDataById({id: item.id, data: item.data}))
            })
            setPrev({...prev, isChanged: false})
        }
    }, [prev])

    function getAllOperators(nodes: Node[]): Node[] {
        return nodes.filter((item: Node) => item.type && operatorNodes.includes(item.type))
    }

    function getIONodes(nodes: Node[], edges: Edge[]): Node[] {
        return nodes.reduce((res: Node[], item: Node) => {
            if (item.type === 'io') {
                const edge: Edge | undefined = edges.find(edge => edge.target === item.id)
                if (edge) {
                    const node: Node | undefined = nodes.find(item => item.id === edge.source)
                    if (node && node.type && operatorNodes.includes(node.type)) res.push(item)
                }
            }
            return res
        }, [])
    }

    function getCalculatedArrayOfNodes(operators: Node[], IONodes: Node[], connectionsOfOperators: Edge[]): Node[] {
        const arrayOfRelations: Node[][] = operators.reduce((res: Node[][], node) => {
            const sortedTargetsIds: string[] = connectionsOfOperators.filter(item => item.source === node.id).sort((a, b) => `${a.sourceHandle}` > `${b.sourceHandle}` ? 1 : -1).map(item => item.target)
            const children: Node[] = sortedTargetsIds.reduce((res: Node[], id) => {
                const node: Node | undefined = [...operators, ...IONodes].find(item => item.id === id)
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
                if (index !== 0 && operatorNodes.includes(node.type as string) && typeof node.data !== "number") {
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