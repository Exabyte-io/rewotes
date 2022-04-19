import {Box} from "@mui/material";
import React from "react";
import {Handle, Position} from 'react-flow-renderer'
import TextArea from "./TextArea";

type Props = {
    id: string,
    data: string
}

const Process: React.FC<Props> = (props: Props) => {
    const {id, data} = props


    return (
        <>
            <Handle type="source" position={Position.Top} id={'a'}/>
            <Handle type="source" position={Position.Bottom} id={'b'}/>
            <Box sx={{
                width: '150px',
                height: '50px',
                border: '1px solid black',
                display: 'flex',
                justifyContent: 'center',
                alignItems: 'center',
                backgroundColor: 'white'
            }}>
                <TextArea id={id} data={data}/>
            </Box>
        </>
    );
}

export default Process