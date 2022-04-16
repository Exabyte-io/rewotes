import {Box} from "@mui/material";
import React from "react";
import {Node, Edge} from 'react-flow-renderer'
import {useSelector} from "react-redux";
import {ReducersType} from "../redux/store";


const SideBar: React.FC = () => {
    const nodes: Node[] = useSelector((state: ReducersType) => state.nodes)
    const edges: Edge[] = useSelector((state: ReducersType) => state.edges)


    function getStringData(arr: Node[] | Edge[]): JSX.Element[] {
        return arr.reduce((res: JSX.Element[], item: Node | Edge) => {
            res.push(
                <li key={item.id} style={{wordBreak: 'break-all'}}>
                    {JSON.stringify(item)}
                </li>
            )
            return res
        }, [])
    }

    return (
        <Box>
            <Box sx={{width: '100%', textAlign: 'start'}}>
                <Box sx={{textAlign: 'center', fontSize: '20px', fontWeight: 'bold'}}>
                    <span>Nodes</span>
                </Box>
                <ol>
                    {getStringData(nodes)}
                </ol>
            </Box>
            <Box sx={{width: '100%', textAlign: 'start'}}>
                <Box sx={{textAlign: 'center', fontSize: '20px', fontWeight: 'bold'}}>
                    <span>Edges</span>
                </Box>
                <ol>
                    {getStringData(edges)}
                </ol>
            </Box>
        </Box>
    )
}

export default SideBar