import React, { useEffect } from 'react';
import hljs from 'highlight.js';

const JSONViewer = ({ children, nodes, edges, isDarkMode }) => {
    // Darkmode style toggling
    const jsonViewerStyle = {
        backgroundColor: isDarkMode
            ? 'rgba(30, 30, 30, 1)'
            : 'rgba(255, 255, 255, 1)',
        color: isDarkMode ? 'rgba(255, 255, 255, 1)' : 'rgba(0, 0, 0, 1)',
    };

    useEffect(() => {
        hljs.highlightAll();
    }, [nodes, edges]);

    return (
        <div
            className='json-viewer'
            style={{ ...jsonViewerStyle, height: '100vh' }}
        >
            {children}
            <div className='json-contents'>
                <pre>
                    nodes: [
                    {nodes.map((node) => {
                        return (
                            <code key={node.id}>
                                {JSON.stringify(
                                    node,
                                    ['id', 'data', 'value', 'type'],
                                    4
                                )}
                            </code>
                        );
                    })}
                    ]
                    <br /><br />
                    edges: [
                    {edges.map((edge) => {
                        return (
                            <code key={edge.id}>
                                {JSON.stringify(
                                    edge,
                                    [
                                        'id',
                                        'source',
                                        'sourceHandle',
                                        'target',
                                        'targetHandle',
                                    ],
                                    4
                                )}
                            </code>
                        );
                    })}
                    ]
                </pre>
            </div>
        </div>
    );
};

export default JSONViewer;
