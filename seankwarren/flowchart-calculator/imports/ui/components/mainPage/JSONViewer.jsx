import React, { useEffect } from 'react';
import hljs from 'highlight.js';

const JSONViewer = ({ 
    children, 
    nodes, 
    edges, 
    flows, 
    loadFlow, 
    isDarkMode, 
    onSave, 
    flowName,
    setFlowName,
}) => {
    // Darkmode style toggling
    const jsonViewerStyle = {
        backgroundColor: isDarkMode
            ? 'rgba(30, 30, 30, 1)'
            : 'rgba(255, 255, 255, 1)',
        color: isDarkMode ? 'rgba(255, 255, 255, 1)' : 'rgba(0, 0, 0, 1)',
    };

    const handleFlowNameChange = (e) => {
        setFlowName(e.target.value);
      };

    useEffect(() => {
        hljs.highlightAll();
    }, [nodes, edges]);

    // TODO: move save/load controls to seperate component
    return (
        <div
            className='json-viewer'
            style={{ ...jsonViewerStyle, height: '100vh' }}
        >
            {children}
            <div className="save-control-panel">
                <select
                    onChange={(e) => loadFlow(flows[e.target.selectedIndex - 1])}
                    style={{ margin: "1rem" }}
                >
                    <option>Select a flow</option>
                    {
                        flows.map((flow) => (
                            <option key={flow._id} value={flow._id}>
                            {flow.name ? `${flow.name}` : flow._id}
                            </option>
                        ))
                    }
                </select>
                <button onClick={onSave}>Save Flow</button>
                <input
                    type="text"
                    placeholder="Flow Name"
                    value={flowName}
                    onChange={handleFlowNameChange}
                />
            </div>
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
