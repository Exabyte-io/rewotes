// Helper function to perform an operation
const performOperation = (op, val1, val2) => {
    if (!val2) {
        switch (op) {
            case 'sin':
                return Math.sin(val1);
            case 'cos':
                return Math.cos(val1);
            case 'tan':
                return Math.tan(val1);
            case 'exp':
                return Math.exp(val1);
        }
    }

    switch (op) {
        case 'add':
            return val1 + val2;
        case 'subtract':
            return val1 - val2;
        case 'multiply':
            return val1 * val2;
        case 'divide':
            return val1 / val2;
        case 'exponent':
            return val1 ** val2;
        case 'greater':
            return val1 > val2;
        case 'less':
            return val1 < val2;
        case 'greaterEqual':
            return val1 >= val2;
        case 'lessEqual':
            return val1 <= val2;
        case 'equal':
            return val1 === val2;
        case 'notEqual':
            return val1 !== val2;
        default:
            return 0;
    }
};

export default performOperation;