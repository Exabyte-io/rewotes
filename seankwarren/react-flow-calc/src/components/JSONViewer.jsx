import React from 'react'

const JSONViewer = ({ children, nodes, edges, isDarkMode }) => {
    
    const jsonViewerStyle = {
        backgroundColor: isDarkMode ? 'rgba(30, 30, 30, 1)' : 'rgba(255, 255, 255, 1)',
        color: isDarkMode ? 'rgba(255, 255, 255, 1)' : 'rgba(0, 0, 0, 1)',
    };

    return (
        <div className="json-viewer" style={{  ...jsonViewerStyle, height: '100vh' }}>
            { children }
            <div className="json-contents">
                nodes: [
                {nodes.map((node) => {
                    return <p key={node.id}>{JSON.stringify(node, ['id', 'data', 'value', 'type'], 4)}</p>
                })}]
                <br />
                edges: [
                {edges.map((edge) => {
                    return <p key={edge.id}>{JSON.stringify(edge, ['id', 'source', 'sourceHandle', 'target', 'targetHandle'], 4)}</p>
                })}]
            </div>
        </div>
    )
}

export default JSONViewer