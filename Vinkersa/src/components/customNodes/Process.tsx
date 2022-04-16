import {Box} from "@mui/material";
import React, {useState} from "react";
import {Handle, Position} from 'react-flow-renderer'

const Process: React.FC = () => {
    const [text, setText] = useState<string>('')

    return (
        <>
            <Handle type="target" position={Position.Top}/>
            <Box sx={{
                minWidth: 150,
                minHeight: 50,
                border: '1px solid black',
                display: 'flex',
                justifyContent: 'center',
                alignItems: 'center',
                backgroundColor: 'white'
            }}>
                <textarea style={{resize: 'none', outline: 'none', border: 'none'}} value={text}
                          onChange={(e) => setText(e.target.value)}/>
            </Box>
            <Handle type="source" position={Position.Bottom} id="a"/>
        </>
    );
}

export default Process