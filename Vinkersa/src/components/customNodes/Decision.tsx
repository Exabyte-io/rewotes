import {Box} from "@mui/material";
import React from "react";
import {Handle, Position} from 'react-flow-renderer'
import TextArea from "./TextArea";

const Decision: React.FC = () => {
    return (
        <>
            <Handle type="target" position={Position.Top}/>
            <Box sx={{
                width: '142px',
                height: '142px',
                display: 'flex',
                justifyContent: 'center',
                alignItems: 'center',
            }}>
                <Box sx={{
                    width: '100px',
                    height: '100px',
                    border: '1px solid black',
                    transform: 'rotate(45deg)',
                    backgroundColor: 'white',
                    display: 'flex',
                    justifyContent: 'center',
                    alignItems: 'center',
                }}>
                    <Box sx={{
                        transform: 'rotate(-45deg)',
                        width: '90px'
                    }}>
                        <TextArea/>
                    </Box>
                </Box>
            </Box>
            <Handle type="source" position={Position.Left} id="a"/>
            <Handle type="source" position={Position.Right} id="b"/>
        </>
    );
}

export default Decision