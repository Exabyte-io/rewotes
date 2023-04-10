import React from 'react';

const DraggableButton = ({ label, nodeType, onDragStart }) => {
    return (
        <button
            className={`node button ${nodeType}`}
            draggable='true'
            onDragStart={(e) => onDragStart(e, nodeType)}
        >
            <span>{label}</span>
        </button>
    );
};

export default DraggableButton;