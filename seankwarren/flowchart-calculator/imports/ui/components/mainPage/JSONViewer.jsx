import React, { useEffect } from 'react';
import hljs from 'highlight.js/lib/core';
import json from 'highlight.js/lib/languages/json';

const JSONViewer = ({
    children,
    nodes,
    edges,
    flows,
    loadFlowchart,
    isDarkMode,
    onSave,
    flowName,
    updateFlowName,
}) => {
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

    useEffect(() => {
        hljs.registerLanguage('json', json);
    }, [])
    

    // TODO: move save/load controls to seperate component
    return (
        <div
            className='json-viewer'
            style={{ ...jsonViewerStyle, height: '100vh' }}
        >
            {children}
            <div className='save-control-panel'>
                <select
                    className='flow-dropdown'
                    onChange={(e) =>
                        loadFlowchart(flows[e.target.selectedIndex - 1])
                    }
                >
                    <option>Select a flow</option>
                    {flows.map((flow) => (
                        <option key={flow._id} value={flow._id}>
                            {flow.name ? `${flow.name}` : flow._id}
                        </option>
                    ))}
                </select>
                <button className='saveflow button' onClick={onSave}>
                    Save Flow
                </button>
                <input
                    className='flowname text-input'
                    type='text'
                    placeholder='Flow Name'
                    value={flowName}
                    onChange={updateFlowName}
                />
            </div>
            <div className='json-contents'>
                <pre>
                    nodes: [
                    {nodes.map((node) => {
                        return (
                            <code className="json" key={node.id}>
                                {JSON.stringify(
                                    node,
                                    ['id', 'data', 'value', 'type'],
                                    4
                                )}
                            </code>
                        );
                    })}
                    ]
                    <br />
                    <br />
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
