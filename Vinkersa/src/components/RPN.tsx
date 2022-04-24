import React, {useEffect, useState} from "react";
import {Edge, Node} from "react-flow-renderer";
import {ReducersType} from "../redux/store";
import {useDispatch, useSelector} from "react-redux";
import {updateNodes} from "../redux/actions";
import {operatorNodes} from "./customNodes/CustomNodeTypes";

const RPN: React.FC = () => {
    const dispatch = useDispatch()

    const nodes: Node[] = useSelector((state: ReducersType) => state.nodes)
    const edges: Edge[] = useSelector((state: ReducersType) => state.edges)

    const [operators, setOperators] = useState<Node[]>([])
    const [IONodes, setIONodes] = useState<Node[]>([])

    const [prev, setPrev] = useState<{state: string, isChanged: boolean}>({state: '', isChanged: false})

    useEffect(() => {
        setOperators(getOperators(nodes, edges))
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
            dispatch(updateNodes(getCalculatedArrayOfNodes(operators, IONodes, edges)))
            setPrev({...prev, isChanged: false})
        }
    }, [prev])

    function getOperators(nodes: Node[], edges: Edge[]): Node[] {
        return nodes.reduce((res: Node[], item: Node) => {
            if (item && item.type && operatorNodes.includes(item.type)) {
                const edge: Edge[] = edges.filter(edge => edge.source === item.id)
                if (edge.length > 0) res.push(item)
            }
            return res
        }, [])
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
                    const isChildOperatorHasAChildren = arr.reduce((res: boolean, item: Node[]) => {
                        if (item[0].id === node.id && item.length > 1) res = true
                        return res
                    }, false)
                    if (isChildOperatorHasAChildren) {
                        isCalculated = false
                        recall = true
                    }
                }
            })
            if (isCalculated) nodes[0].data = calculateValues(nodes)
        })

        return recall ? calculateNodes(arr) : arr.map(item => item[0])
    }

    function calculateValues(nodes: Node[]): number {
        const typeOfOperator = nodes[0].type
        return nodes.reduce((res: number, item: Node, index: number) => {
            const data = item.data === '' ? 0 : + item.data
            if (index === 1) res = data
            if (index !== 0 && index !== 1) {
                switch (typeOfOperator) {
                    case 'plus': return res += data
                    case 'minus': return res -= data
                    case 'divide': return res /= data
                    case 'multiply': return res *= data
                }
            }
            return res
        }, 0)
    }

    //TODO add view for rpn
    return <></>
}

export default RPN