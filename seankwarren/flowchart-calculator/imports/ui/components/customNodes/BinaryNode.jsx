import React from 'react';
import { Handle } from 'reactflow';
import { binaryOperations } from '../../utils/operationDef';

const BinaryNode = ({ id, data }) => {
    const handleSelectChange = (e) => {
        data.onChange(e.target.value);
    };

    return (
        <div className='node binary round' data-testid='binary-node'>
            <Handle
                className='handle input top'
                type='target'
                position='left'
                id={`${id}-left-top`}
                style={{
                    top: '30%',
                }}
            />
            <Handle
                className='handle input bottom'
                type='target'
                position='left'
                id={`${id}-left-bottom`}
                style={{
                    top: '70%',
                }}
            />
            <select
                className='node-dropdown'
                defaultValue={data.value ? data.value : '+'}
                onChange={handleSelectChange}
            >
                {binaryOperations.map((op) => (
                    <option key={op.value} value={op.value}>
                        {op.content}
                    </option>
                ))}
            </select>
            <Handle
                className='handle output'
                type='source'
                position='right'
                id={`${id}-right`}
            />
        </div>
    );
};

export default BinaryNode;
