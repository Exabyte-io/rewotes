import React from 'react'
import { Handle } from 'reactflow';

const UnaryNode = ({ id, data }) => {

    const handleSelectChange = (e) => {
        data.onChange(e.target.value);
    };

    return (
        <div className='node unary round'>
            <Handle
                className='handle input'
                type='target'
                position='left'
                id={`${id}-left`}
            />
            <select onChange={handleSelectChange}>
                <option value='sin'>sin</option>
                <option value='cos'>cos</option>
                <option value='tan'>tan</option>
                <option value='exp'>e^x</option>
            </select>
            <Handle
                className='handle output'
                type='source'
                position='right'
                id={`${id}-right`}
            />
        </div>
    );
}

export default UnaryNode