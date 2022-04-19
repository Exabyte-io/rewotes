import {Box} from "@mui/material";
import React from "react";
import {Handle, Position} from 'react-flow-renderer'
import TextArea from "./TextArea";

type Props = {
    id: string
}

const Decision: React.FC<Props> = (props: Props) => {
    const {id} = props

    return (
        <>
            <Handle type="target" position={Position.Top}/>
            <Handle type="source" position={Position.Left} id="a"/>
            <Handle type="source" position={Position.Right} id="b"/>
            <Handle type="source" position={Position.Bottom} id="c"/>
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
                        <TextArea id={id}/>
                    </Box>
                </Box>
            </Box>
        </>
    );
}

export default Decision