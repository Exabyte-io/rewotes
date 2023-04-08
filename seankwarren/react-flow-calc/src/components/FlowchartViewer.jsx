import React from 'react'
import ReactFlow, {
    Controls,
    Background,
    MiniMap,
} from 'reactflow'

const FlowchartViewer = () => {
    return (
        <ReactFlow>
            <MiniMap />
            <Controls />
            <Background />
        </ReactFlow>
    )
}

export default FlowchartViewer