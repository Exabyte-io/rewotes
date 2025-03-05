import React, { useEffect } from 'react';
import hljs from 'highlight.js/lib/core';
import json from 'highlight.js/lib/languages/json';
import { useDarkMode } from '../reusable/DarkModeContext';

const JSONViewer = ({
    children,
    nodes,
    edges,
    flows,
    loadFlowchart,
    onSave,
    flowName,
    updateFlowName,
}) => {

    const { isDarkMode } = useDarkMode();

    useEffect(() => {
        hljs.registerLanguage('json', json);
    }, [])
    
    useEffect(() => {
        hljs.highlightAll();
    }, [nodes, edges]);
    

    // TODO: move save/load controls to seperate component
    return (
        <div
            className={'json-viewer'+ (isDarkMode ? " dark-mode" : "")}
            style={{ height: '100vh' }}
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
