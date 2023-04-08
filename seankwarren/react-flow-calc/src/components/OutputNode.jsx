import React from 'react'
import { Handle } from 'reactflow';

const OutputNode = ({ id, data}) => {
  return (
    <div className='node output' data-testid='output-node'>
    <Handle
        type='target'
        position='left'
        id={`${id}-left`}
        style={{ background: '#ffffff', borderColor: '#000000' }}
        // onConnect={(params) => console.log('handle onConnect', params)}
    />
    <strong>Result: </strong> {String(data.value)}
</div>
  )
}

export default OutputNode