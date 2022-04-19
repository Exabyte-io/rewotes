import {Box} from "@mui/material";
import React from "react";
import {Handle, Position} from 'react-flow-renderer'
import TextArea from "./TextArea";

type Props = {
    id: string
}

const IO: React.FC<Props> = (props: Props) => {
    const {id} = props


    return (
        <>
            <Handle type="target" position={Position.Top}/>
            <Handle type="source" position={Position.Bottom} id="a"/>
            <Box sx={{
                width: '150px',
                height: '50px',
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
                    <TextArea id={id}/>
                </Box>
            </Box>
        </>
    );
}

export default IO