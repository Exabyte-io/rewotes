import {Box} from "@mui/material";
import React from "react";
import {Handle, Position} from 'react-flow-renderer'
import TextArea from "./TextArea";

type Props = {
    id: string,
    data: string
}

const IO: React.FC<Props> = (props: Props) => {
    const {id, data} = props


    return (
        <>
            <Handle type="target" position={Position.Top} id={'a'}/>
            <Box sx={{
                width: '150px',
                height: '50px',
                transform: 'skew(-25deg)',
                border: '1px solid black',
                display: 'flex',
                justifyContent: 'center',
                alignItems: 'center',
                backgroundColor: 'white'
            }}/>
            <Box style={{position: 'absolute', top: '8px', left: '10px', width: '80%'}}>
                <TextArea id={id} data={data}/>
            </Box>

        </>
    );
}

export default IO