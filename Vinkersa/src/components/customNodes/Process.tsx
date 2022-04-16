import {Box} from "@mui/material";
import React from "react";
import {Handle, Position} from 'react-flow-renderer'
import TextArea from "./TextArea";

const Process: React.FC = () => {
    return (
        <>
            <Handle type="target" position={Position.Top}/>
            <Box sx={{
                width: '150px',
                height: '50px',
                border: '1px solid black',
                display: 'flex',
                justifyContent: 'center',
                alignItems: 'center',
                backgroundColor: 'white'
            }}>
                <TextArea/>
            </Box>
            <Handle type="source" position={Position.Bottom} id="a"/>
        </>
    );
}

export default Process