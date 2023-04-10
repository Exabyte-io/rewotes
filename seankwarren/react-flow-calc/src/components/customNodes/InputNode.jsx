import React, { useState } from 'react';
import { Handle } from 'reactflow';

const InputNode = ({ id, data }) => {
    const handleInputChange = (e) => {
        const inputValue = parseFloat(e.target.value);
        data.onChange(inputValue);
    };

    return (
        <div className='node input' data-testid='input-node'>
            <input type='number' defaultValue={0} onChange={handleInputChange} />
            <Handle
                className='handle output'
                type='source'
                position='right'
                id={`${id}-right`}
            />
        </div>
    );
};

export default InputNode;
