import {Box} from "@mui/material";
import React from "react";
import {Handle, Position} from 'react-flow-renderer'
import TextArea from "./TextArea";

type Props = {
    id: string,
    data: string,
    isOperator?: boolean,
    operator?: string
}

const Process: React.FC<Props> = (props: Props) => {
    const {id, data, isOperator, operator} = props


    return (
        <>
            <Handle type="source" position={Position.Top} id={'a'}/>
            {isOperator ?
                <Box sx={{
                    position: 'absolute',
                    bottom: '-4px',
                    display: 'flex',
                    justifyContent: 'space-around',
                    width: '100%'
                }}>
                    <Handle type="source" position={Position.Bottom} style={{position: 'static'}} id={'b'}/>
                    <Handle type="source" position={Position.Bottom} style={{position: 'static'}} id={'c'}/>
                </Box>
                :
                <Handle type="source" position={Position.Bottom} id={'b'}/>
            }
            <Box sx={{
                width: '150px',
                height: '50px',
                border: '1px solid black',
                display: 'flex',
                justifyContent: 'center',
                alignItems: 'center',
                backgroundColor: 'white'
            }}>
                <TextArea id={id} data={data} isReadOnly={isOperator} prefix={operator}/>
            </Box>
        </>
    );
}

export default Process