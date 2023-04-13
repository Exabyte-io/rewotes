import React from 'react';

const DraggableButton = ({ label, nodeType, onDragStart, className }) => {
    return (
        <button
            className={`node button ${nodeType} ${className}`}
            draggable='true'
            onDragStart={(e) => onDragStart(e, nodeType)}
        >
            <span>{label}</span>
        </button>
    );
};

export default DraggableButton;
