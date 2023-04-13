import React from 'react';
import DraggableButton from '../reusable/DraggableButton';

const NodeButtons = ({ handleDragStart, isDarkMode, clearFlowchart, clearFlows }) => {
    // Darkmode style toggling
    const buttonsPanelStyle = {
        backgroundColor: isDarkMode
            ? 'rgba(0, 0, 0, 0.5)'
            : 'rgba(51, 51, 51, 0.3)',
    };

    return (
        <div className='buttons-panel' style={{ ...buttonsPanelStyle, width: '100%' }}>
            <DraggableButton label="in" nodeType="input" onDragStart={handleDragStart} />
            <DraggableButton label="+" nodeType="binary" onDragStart={handleDragStart} className="round"/>
            <DraggableButton label="sin" nodeType="unary" onDragStart={handleDragStart} className="round"/>
            <DraggableButton label="&gt;" nodeType="comparison" onDragStart={handleDragStart} className="round"/>
            <DraggableButton label="out" nodeType="output" onDragStart={handleDragStart} />
            <button className="clear" onClick={clearFlowchart}>Clear</button>
            {/* <button className="empty" onClick={clearFlows}>Empty</button> */}
        </div>
    );
};

export default NodeButtons;
