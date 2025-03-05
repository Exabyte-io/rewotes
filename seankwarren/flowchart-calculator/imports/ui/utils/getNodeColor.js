const getNodeColor = (node) => {
    switch (node.type) {
        case 'inputNode':
            return 'var(--input-color)';
        case 'binaryNode':
            return 'var(--binary-color)';
        case 'unaryNode':
            return 'var(--unary-color)';
        case 'comparisonNode':
            return 'var(--comparison-color)';
        case 'outputNode':
            return 'var(--output-color)';
        default:
            return '#eee';
    }
};

export default getNodeColor;