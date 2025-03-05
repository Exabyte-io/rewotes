import React from 'react';
import { Handle } from 'reactflow';

const OutputNode = ({ id, data }) => {
    return (
        <div className='node output' data-testid='output-node' data-nodeid={id}>
            <Handle
                className='handle input'
                type='target'
                position='left'
                id={`${id}-left`}
            />
            <strong>Result: </strong> {String(data.value)}
        </div>
    );
};

export default OutputNode;
