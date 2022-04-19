import React from "react";
import {Node} from "react-flow-renderer";
import {useDispatch, useSelector} from "react-redux";
import {ReducersType} from "../../redux/store";
import {setNodesData} from "../../redux/actions";

type Props = {
    id: string
}

const TextArea: React.FC<Props> = (props: Props) => {
    const {id} = props
    const nodes: Node[] = useSelector((state: ReducersType) => state.nodes)
    const dispatch = useDispatch()

    function getValue(nodes: Node[], id: string): string {
        const node: Node | undefined = nodes.find(item => item.id === id)
        if (node && node.data) return node.data
        return ''
    }
    console.log('123')
    function changeData(nodes: Node[], id: string, text: string) {
        const index: number = nodes.findIndex(item => item.id === id)
        if (index >= 0 && nodes[index]) {
            nodes[index].data = text
            dispatch(setNodesData([...nodes]))
        }
    }

    return (
        <textarea style={{resize: 'none', outline: 'none', border: 'none', width: '100%'}} value={getValue(nodes, id)}
                  onChange={(e) => changeData(nodes, id, e.target.value)}/>
    )
}

export default TextArea