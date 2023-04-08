import React from 'react'

const NodeButtons = ({ addNode, handleDragStart, isDarkMode }) => {
    const buttonsPanelStyle = {
        backgroundColor: isDarkMode ? 'rgba(0, 0, 0, 0.5)' : 'rgba(51, 51, 51, 0.3)',
    };

    return (
        <div className='buttons-panel' style={{ ...buttonsPanelStyle, width: '100%'}}>
            <button
                className='node button input'
                draggable='true'
                onDragStart={(e) => handleDragStart(e, 'input')}
            >
                in
            </button>
            <button
                className='node button binary round'
                draggable='true'
                onDragStart={(e) => handleDragStart(e, 'binary')}
            >
                +
            </button>
            <button
                className='node button unary round'
                draggable='true'
                onDragStart={(e) => handleDragStart(e, 'unary')}
            >
                sin
            </button>
            <button
                className='node button comparison round'
                draggable='true'
                onDragStart={(e) => handleDragStart(e, 'comparison')}
            >
                &gt;
            </button>
            <button
                className='node button output'
                draggable='true'
                onDragStart={(e) => handleDragStart(e, 'output')}
            >
                out
            </button>
        </div>
    );
};

export default NodeButtons;