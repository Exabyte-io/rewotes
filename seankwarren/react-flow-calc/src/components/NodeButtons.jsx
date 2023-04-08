import React from 'react'

const NodeButtons = ({ addNode, handleDragStart }) => {
    return (
        <div className='buttons-panel'>
            <button
                className='node-button input'
                draggable='true'
                onDragStart={(e) => handleDragStart(e, 'input')}
            >
                in
            </button>
            <button
                className='node-button operation round'
                draggable='true'
                onDragStart={(e) => handleDragStart(e, 'operation')}
            >
                +
            </button>
            <button
                className='node-button comparison round'
                draggable='true'
                onDragStart={(e) => handleDragStart(e, 'comparison')}
            >
                &gt;
            </button>
            <button
                className='node-button output'
                draggable='true'
                onDragStart={(e) => handleDragStart(e, 'output')}
            >
                out
            </button>
        </div>
    );
};

export default NodeButtons;