import React, { useState } from 'react';

const useDraggable = (reactFlowInstance, addNode) => {
    // Drag and drop state
    const [draggedNodeType, setDraggedNodeType] = useState(null);

    const handleDragStart = (e, nodeType) => {
        setDraggedNodeType(nodeType);
        e.dataTransfer.setData('text/plain', nodeType);
    };

    const handleDragOver = (e) => {
        e.preventDefault();
        e.stopPropagation();
    };

    const handleDrop = (e) => {
        e.preventDefault();
        e.stopPropagation();

        if (!draggedNodeType) return;

        const reactFlowBounds = e.currentTarget.getBoundingClientRect();
        const position = reactFlowInstance.project({
            x: e.clientX - reactFlowBounds.left,
            y: e.clientY - reactFlowBounds.top,
        });

        addNode(draggedNodeType, position);
        setDraggedNodeType(null);
    };

    return { handleDragStart, handleDragOver, handleDrop };
};

export default useDraggable;
