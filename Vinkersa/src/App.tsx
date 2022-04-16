import { Box, Button, MenuItem, TextField } from "@mui/material";
import React, {useEffect, useState} from "react";
import FlowChart from "./pages/FlowChart";
import {Node, Edge} from 'react-flow-renderer'
import SideBar from "./components/SideBar";
import { useDispatch } from "react-redux";
import {setEdgesData, setNodesData} from "./redux/actions";
import {NodeNameTypes} from "./components/customNodes/CustomNodeTypes";
import {store} from "./redux/store";
import {nanoid} from 'nanoid'

const initialNodes: Node[] = [
    {
        id: '1',
        type: 'terminal',
        data: null,
        position: { x: 250, y: 25 },
    },

    {
        id: '2',
        type: 'io',
        data: null,
        position: { x: 100, y: 125 },
    },
    {
        id: '3',
        type: 'process',
        data: null,
        position: { x: 250, y: 250 },
    },
    {
        id: '4',
        type: 'decision',
        data: null,
        position: { x: 250, y: 350 },
    }
];

const initialEdges: Edge[] = [
    { id: 'e1-2', source: '1', target: '2'},
    { id: 'e2-3', source: '2', target: '3'},
    { id: 'e3-4', source: '3', target: '4'},
];

const App: React.FC = () => {
    const dispatch = useDispatch()
    const [nameOfNode, setNameOfNode] = useState<NodeNameTypes | ''>('')

    useEffect(() => {
        dispatch(setNodesData(initialNodes))
        dispatch(setEdgesData(initialEdges))
    }, [])

    function clickHandler(nameOfNode: NodeNameTypes) {
        setNameOfNode('')
        dispatch(setNodesData([...store.getState().nodes, {
            id: nanoid(8),
            type: nameOfNode,
            data: null,
            position: {x: 300, y: 300}
        }]))
    }


  return (
      <Box sx={{width: '100%', height: '100%', display: 'flex'}}>
          <Box sx={{width: '80%', height: '100%'}}>
              <FlowChart/>
          </Box>
          <Box sx={{width: '20%', height: '100%'}}>
                <Box sx={{width: '100%', display: 'flex', justifyContent: 'space-around', alignItems: 'center', paddingTop: '20px', paddingBottom: '20px'}}>
                    <TextField
                        sx={{width: '50%'}}
                        select
                        label="Node"
                        value={nameOfNode}
                        onChange={(e) => setNameOfNode(e.target.value as NodeNameTypes)}
                    >
                        <MenuItem value={'terminal'}>Terminal</MenuItem>
                        <MenuItem value={'io'}>I/O</MenuItem>
                        <MenuItem value={'process'}>Process</MenuItem>
                        <MenuItem value={'decision'}>Decision</MenuItem>
                    </TextField>
                    <Button variant='contained' disabled={nameOfNode === ''} onClick={() => clickHandler(nameOfNode as NodeNameTypes)}>Add block</Button>
                </Box>
              <SideBar/>
          </Box>
      </Box>
  )
}

export default App
