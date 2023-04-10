import React from 'react';
import { Handle } from 'reactflow';

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
            <select onChange={handleSelectChange}>
                <option value='greater'>{'>'}</option>
                <option value='less'>{'<'}</option>
                <option value='greaterEqual'>{'>='}</option>
                <option value='lessEqual'>{'<='}</option>
                <option value='equal'>{'=='}</option>
                <option value='notEqual'>{'!='}</option>
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
