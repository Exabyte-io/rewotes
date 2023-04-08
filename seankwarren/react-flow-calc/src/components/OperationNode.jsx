import React from 'react'
import { Handle } from 'reactflow';

const OperationNode = ({ id, data }) => {

    const handleSelectChange = (e) => {
        data.onChange(e.target.value);
    };

    return (
        <div className='node operation' data-testid='operation-node'>
            <Handle
                type='target'
                position='left'
                id={`${id}-left-top`}
                style={{
                    background: '#ffffff',
                    borderColor: '#000000',
                    top: '30%',
                }}
                // onConnect={(params) => console.log('handle onConnect', params)}
            />
            <Handle
                type='target'
                position='left'
                id={`${id}-left-bottom`}
                style={{
                    background: '#ffffff',
                    borderColor: '#000000',
                    top: '70%',
                }}
                // onConnect={(params) => console.log('handle onConnect', params)}
            />
            <select onChange={handleSelectChange}>
                <option value='add'>+</option>
                <option value='subtract'>-</option>
                <option value='multiply'>*</option>
                <option value='divide'>/</option>
            </select>
            <Handle
                type='source'
                position='right'
                id={`${id}-right`}
                style={{ background: '#000000' }}
                // onConnect={(params) => console.log('handle onConnect', params)}
            />
        </div>
    );
}

export default OperationNode