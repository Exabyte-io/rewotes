import {Box} from "@mui/material";
import React, {useState} from "react";
import {Handle, Position} from 'react-flow-renderer'

const IO: React.FC = () => {
    const [text, setText] = useState<string>('')

    return (
        <>
            <Handle type="target" position={Position.Top}/>
            <Box sx={{
                minWidth: 150,
                minHeight: 50,
                transform: 'skew(-25deg)',
                border: '1px solid black',
                display: 'flex',
                justifyContent: 'center',
                alignItems: 'center',
                backgroundColor: 'white'
            }}>
                <Box sx={{
                    transform: 'skew(25deg)',
                    paddingLeft: 2,
                    paddingRight: 2,
                }}>
                    <textarea style={{resize: 'none', outline: 'none', border: 'none'}} value={text}
                              onChange={(e) => setText(e.target.value)}/>
                </Box>
            </Box>
            <Handle type="source" position={Position.Bottom} id="a"/>
        </>
    );
}

export default IO