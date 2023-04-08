import React from 'react'

const NodeButtons = ({ addNode, handleDragStart }) => {
    return (
        <div className='node-buttons'>
            <button
                className='node-button input-button'
                draggable='true'
                onDragStart={(e) => handleDragStart(e, 'input')}
            >
                Add Input
            </button>
            <button
                className='node-button operation-button'
                draggable='true'
                onDragStart={(e) => handleDragStart(e, 'operation')}
            >
                Add Operation
            </button>
            <button
                className='node-button output-button'
                draggable='true'
                onDragStart={(e) => handleDragStart(e, 'output')}
            >
                Add Output
            </button>
            <button
                className='node-button comparison-button'
                draggable='true'
                onDragStart={(e) => handleDragStart(e, 'comparison')}
            >
                Add Comparison
            </button>
        </div>
    );
};

export default NodeButtons;