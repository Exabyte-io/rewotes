import React from 'react';
import {
    InputNode,
    BinaryNode,
    UnaryNode,
    OutputNode,
    ComparisonNode,
} from './index';

const nodeTypesConfig = {
    inputNode: (props) => <InputNode {...props} />,
    binaryNode: (props) => <BinaryNode {...props} />,
    unaryNode: (props) => <UnaryNode {...props} />,
    comparisonNode: (props) => <ComparisonNode {...props} />,
    outputNode: (props) => <OutputNode {...props} />,
};

export default nodeTypesConfig;
