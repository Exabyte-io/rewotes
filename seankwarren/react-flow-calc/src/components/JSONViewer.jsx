import React from 'react'

const JSONViewer = ({ nodes, edges }) => {
    return (
        <div className="json-viewer" style={{ height: '100vh' }}>
            nodes:
            {nodes.map((node) => {
                return <p>{JSON.stringify(node, ['id', 'data', 'value', 'type'], 4)}</p>
            })}
            <br />
            edges:
            {edges.map((edge) => {
                return <p>{JSON.stringify(edge, ['id', 'source', 'sourceHandle', 'target', 'targetHandle'], 4)}</p>
            })}
        </div>
    )
}

export default JSONViewer