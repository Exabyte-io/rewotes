import React from 'react';
import { Handle } from 'reactflow';
import { comparisonOperations } from '../../utils/operationLabels'

export const ComparisonNode = ({ id, data }) => {
    const handleSelectChange = (e) => {
        data.onChange(e.target.value,);
    };

    return (
        <div className='node comparison round' data-testid='comparison-node'>
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
            <select className="node-dropdown" defaultValue={data.value ? data.value : ">"} onChange={handleSelectChange}>
                {comparisonOperations.map((op) => <option key={op.value} value={op.value}>{op.content}</option>)}
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

export default ComparisonNode;
