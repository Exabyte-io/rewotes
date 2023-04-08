import React, { useEffect, useState } from 'react'
import { Handle } from 'reactflow';

const InputNode = ({ id, data }) => {
    const [val, setVal] = useState(data);

    const handleInputChange = (e) => {
        const inputValue = parseFloat(e.target.value);
        data.onChange({ id, data: { ...data, value: inputValue } });
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
                // onConnect={(params) => console.log('handle onConnect', params)}
            />
        </div>
    )
}

export default InputNode