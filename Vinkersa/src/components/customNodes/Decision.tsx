import {Box} from "@mui/material";
import React, {useState} from "react";
import {Handle, Position} from 'react-flow-renderer'

const Decision: React.FC = () => {
    const [text, setText] = useState<string>('')

    return (
        <>
            <Handle type="target" position={Position.Top}/>
            <Box sx={{
                width: '108px',
                height: '108px',
                display: 'flex',
                justifyContent: 'center',
                alignItems: 'center',
            }}>
                <Box sx={{
                    width: '75px',
                    height: '75px',
                    border: '1px solid black',
                    transform: 'rotate(45deg)',
                }}>

                </Box>
            </Box>
            <Handle type="source" position={Position.Left} id="a"/>
            <Handle type="source" position={Position.Right} id="b"/>
        </>
    );
}

export default Decision