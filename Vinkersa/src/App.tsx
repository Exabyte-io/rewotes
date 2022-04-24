import {Box, Button, MenuItem, TextField} from "@mui/material";
import React, {useEffect, useState} from "react";
import FlowChart from "./pages/FlowChart";
import {Node, Edge} from 'react-flow-renderer'
import SideBar from "./components/SideBar";
import {useDispatch} from "react-redux";
import {setEdgesData, setNodesData} from "./redux/actions";
import {baseNodes, NodeNameTypes, operatorNodes} from "./components/customNodes/CustomNodeTypes";
import {store} from "./redux/store";
import {nanoid} from 'nanoid'
import RPN from "./components/RPN";


const App: React.FC = () => {
    const dispatch = useDispatch()
    const [nameOfNode, setNameOfNode] = useState<NodeNameTypes | ''>('')

    useEffect(() => {
        const nodes: Node[] = localStorage.getItem('nodes') ? JSON.parse(localStorage.getItem('nodes')) : []
        const edges: Edge[] = localStorage.getItem('edges') ? JSON.parse(localStorage.getItem('edges')) : []
        if (nodes.length === 0) {
            localStorage.setItem('nodes', '')
            localStorage.setItem('edges', '')
        } else if (nodes.length !== 0) {
            dispatch(setNodesData(nodes))
            dispatch(setEdgesData(edges))
        }
    }, [])

    function clickHandler(nameOfNode: NodeNameTypes) {
        setNameOfNode('')
        dispatch(setNodesData([...store.getState().nodes, {
            id: nanoid(8),
            type: nameOfNode,
            data: '',
            position: {x: 300, y: 300}
        }]))
    }


    return (
        <Box sx={{width: '100%', height: '100%', display: 'flex'}}>
            <Box sx={{width: '80%', height: '100%'}}>
                <FlowChart/>
            </Box>
            <Box sx={{width: '20%', height: '100%', overflow: 'auto'}}>
                <Box sx={{
                    width: '100%',
                    display: 'flex',
                    justifyContent: 'space-around',
                    alignItems: 'center',
                    paddingTop: '20px',
                    paddingBottom: '20px'
                }}>
                    <TextField
                        sx={{width: '50%'}}
                        select
                        label="Node"
                        value={nameOfNode}
                        onChange={(e) => setNameOfNode(e.target.value as NodeNameTypes)}
                    >
                        {[...baseNodes, ...operatorNodes].map((item, index) =>
                            <MenuItem key={index} value={item}><span style={{textTransform: 'capitalize'}}>{item}</span></MenuItem>
                        )}
                    </TextField>
                    <Button variant='contained' disabled={nameOfNode === ''}
                            onClick={() => clickHandler(nameOfNode as NodeNameTypes)}>Add block</Button>
                </Box>
                <SideBar/>
                <RPN/>
            </Box>
        </Box>
    )
}

export default App
