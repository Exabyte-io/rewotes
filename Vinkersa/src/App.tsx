import { Box } from "@mui/material";
import React, {useEffect} from "react";
import FlowChart from "./pages/FlowChart";
import {Node, Edge} from 'react-flow-renderer'
import SideBar from "./components/SideBar";
import { useDispatch, useSelector } from "react-redux";
import {ReducersType} from "./redux/store";
import {setEdgesData, setNodesData} from "./redux/actions";

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
    const nodes: Node[] = useSelector((state: ReducersType) => state.nodes)
    const edges: Edge[] = useSelector((state: ReducersType) => state.edges)

    useEffect(() => {
        dispatch(setNodesData(initialNodes))
        dispatch(setEdgesData(initialEdges))
    }, [])


  return (
      <Box sx={{width: '100%', height: '100%', display: 'flex'}}>
          <Box sx={{width: '80%', height: '100%'}}>
              <FlowChart/>
          </Box>
          <Box sx={{width: '20%', height: '100%'}}>
              <SideBar/>
          </Box>
      </Box>
  )
}

export default App
