import {Box} from "@mui/material";
import React from "react";
import {Handle, Position} from 'react-flow-renderer'
import TextArea from "./TextArea";

const Terminal: React.FC = () => {
    return (
        <>
            <Handle type="target" position={Position.Top}/>
            <Box sx={{
                width: '150px',
                height: '50px',
                border: '1px solid black',
                borderRadius: 50,
                paddingLeft: 2,
                paddingRight: 2,
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

export default Terminal