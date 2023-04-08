import React from 'react'

const JSONViewer = ({ nodes, edges }) => {
    return (
        <div style={{ height: '100vh', width: '30vw' }}>
            nodes: {JSON.stringify(nodes, ['id', 'data', 'value', 'type'], 4)}
            edges: {JSON.stringify(edges, ['id', 'source', 'sourceHandle', 'target', 'targetHandle'], 4)}
        </div>
    )
}

export default JSONViewer