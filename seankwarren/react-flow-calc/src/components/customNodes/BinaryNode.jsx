import React from 'react'
import { Handle } from 'reactflow';

const BinaryNode = ({ id, data }) => {

    const handleSelectChange = (e) => {
        data.onChange(e.target.value);
    };

    return (
        <div className='node binary round'>
            <Handle
                className='handle input'
                type='target'
                position='left'
                id={`${id}-left-top`}
                style={{
                    top: '30%',
                }}
            />
            <Handle
                className='handle input'
                type='target'
                position='left'
                id={`${id}-left-bottom`}
                style={{
                    top: '70%',
                }}
            />
            <select onChange={handleSelectChange}>
                <option value='add'>+</option>
                <option value='subtract'>-</option>
                <option value='multiply'>*</option>
                <option value='divide'>/</option>
                <option value='exponent'>^</option>
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

export default BinaryNode