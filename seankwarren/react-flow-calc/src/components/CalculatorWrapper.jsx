import React from 'react'
import FlowchartViewer from './FlowchartViewer'
import JSONViewer from './JSONViewer'
import 'reactflow/dist/style.css';

const CalculatorWrapper = () => {
    return (
        <div style={{ display: 'flex' }}>
            <div style={{ height: '100vh', width: '70vw' }}>
                <FlowchartViewer/>
            </div>
            <JSONViewer />
        </div>
    )
}

export default CalculatorWrapper