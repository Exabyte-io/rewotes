import { Box } from "@mui/material";
import React from "react";
import {Node, Edge} from 'react-flow-renderer'

type Props = {
    nodes: Node[],
    edges: Edge[]
}

const SideBar: React.FC<Props> = (props: Props) => {
    const {nodes, edges} = props

    function getStringData(arr: Node[] | Edge[]): JSX.Element[] {
        return arr.reduce((res: JSX.Element[], item: Node | Edge) => {
            res.push(
                <li key={item.id}>
                    {JSON.stringify(item)}
                </li>
            )
            return res
        }, [])
    }

    return (
        <Box>
            <Box sx={{width: '100%', textAlign: 'center'}}>
                <span>Nodes</span>
            </Box>
            <ol>
                {getStringData(nodes)}
            </ol>
            <Box sx={{width: '100%', textAlign: 'center'}}>
                <span>Edges</span>
            </Box>
            <ol>
                {getStringData(edges)}
            </ol>
        </Box>
    )
}

export default SideBar