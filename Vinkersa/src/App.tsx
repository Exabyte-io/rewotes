import { Box } from "@mui/material";
import React from "react";
import FlowChart from "./pages/FlowChart";
import {Node, Edge} from 'react-flow-renderer'
import SideBar from "./components/SideBar";

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


  return (
      <Box sx={{width: '100%', height: '100%', display: 'flex'}}>
          <Box sx={{width: '80%', height: '100%'}}>
              <FlowChart initialNodes={initialNodes} initialEdges={initialEdges}/>
          </Box>
          <Box sx={{width: '20%', height: '100%'}}>
              <SideBar nodes={initialNodes} edges={initialEdges}/>
          </Box>
      </Box>
  )
}

export default App
