import React from 'react';
import InputNode from './InputNode';
import BinaryNode from './BinaryNode';
import UnaryNode from './UnaryNode';
import OutputNode from './OutputNode';
import ComparisonNode from './ComparisonNode';

const nodeTypesConfig = {
    inputNode: (props) => <InputNode {...props} />,
    binaryNode: (props) => <BinaryNode {...props} />,
    unaryNode: (props) => <UnaryNode {...props} />,
    comparisonNode: (props) => <ComparisonNode {...props} />,
    outputNode: (props) => <OutputNode {...props} />,
};

export default nodeTypesConfig;