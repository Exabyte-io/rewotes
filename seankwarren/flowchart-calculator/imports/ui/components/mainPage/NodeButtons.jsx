import React from 'react';
import DraggableButton from '../reusable/DraggableButton';
import { useDarkMode } from '../reusable/DarkModeContext';

const NodeButtons = ({
    handleDragStart,
    clearFlowchart,
    clearFlows,
}) => {
    const { isDarkMode } = useDarkMode();

    return (
        <div
            className={'buttons-panel'+ (isDarkMode ? " dark-mode" : "")}
            style={{ width: '100%' }}
        >
            <DraggableButton
                label='in'
                nodeType='input'
                onDragStart={handleDragStart}
            />
            <DraggableButton
                label='+'
                nodeType='binary'
                onDragStart={handleDragStart}
                className='round'
            />
            <DraggableButton
                label='sin'
                nodeType='unary'
                onDragStart={handleDragStart}
                className='round'
            />
            <DraggableButton
                label='&gt;'
                nodeType='comparison'
                onDragStart={handleDragStart}
                className='round'
            />
            <DraggableButton
                label='out'
                nodeType='output'
                onDragStart={handleDragStart}
            />
            <button className='clear' onClick={clearFlowchart}>
                Clear
            </button>
            {/* <button className="empty" onClick={clearFlows}>Empty</button> */}
        </div>
    );
};

export default NodeButtons;
