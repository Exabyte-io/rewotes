import React, { useState } from 'react';
import { Handle } from 'reactflow';

const InputNode = ({ id, data }) => {
    const [val, setVal] = useState(data);

    const handleInputChange = (e) => {
        const inputValue = parseFloat(e.target.value)
        data.onChange(inputValue);
        // updateOutputNodes();
        setVal(inputValue);
    };

    return (
        <div className='node input' data-testid='input-node'>
            <input type='number' onChange={handleInputChange} />
            <Handle
                type='source'
                position='right'
                id={`${id}-right`}
                style={{ background: '#000000' }}
            />
        </div>
    );
};

export default InputNode;
