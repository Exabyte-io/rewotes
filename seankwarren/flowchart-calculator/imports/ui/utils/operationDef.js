export const comparisonOperations = [
    {
        value: 'greater',
        content: '>',
        expression: (val1, val2) => val1 > val2,
    },
    {
        value: 'less',
        content: '<',
        expression: (val1, val2) => val1 < val2,
    },
    {
        value: 'greaterEqual',
        content: '>=',
        expression: (val1, val2) => val1 >= val2
    },
    {
        value: 'lessEqual',
        content: '<=',
        expression: (val1, val2) => val1 <= val2
    },
    {
        value: 'equal',
        content: '==',
        expression: (val1, val2) => val1 === val2
    },
    {
        value: 'notEqual',
        content: '!=',
        expression: (val1, val2) => val1 !== val2
    }
]

export const binaryOperations = [
    {
        value: 'add',
        content: '+',
        expression: (val1, val2) => val1 + val2
    },
    {
        value: 'subtract',
        content: '-',
        expression: (val1, val2) => val1 - val2
    },
    {
        value: 'multiply',
        content: '*',
        expression: (val1, val2) => val1 * val2
    },
    {
        value: 'divide',
        content: '/',
        expression: (val1, val2) => val1 / val2
    },
    {
        value: 'exponent',
        content: '^',
        expression: (val1, val2) => val1 ** val2
    },
    {
        value: 'custom',
        content: '(a+b)^2',
        expression: (val1, val2) => (val1 + val2) ** 2
    },
]

export const unaryOperations = [
    {
        value: 'sin',
        content: 'sin',
        expression: (val1, val2) => Math.sin(val1)
    },
    {
        value: 'cos',
        content: 'cos',
        expression: (val1, val2) => Math.cos(val1)
    },
    {
        value: 'tan',
        content: 'tan',
        expression: (val1, val2) => Math.tan(val1)
    },
    {
        value: 'exp',
        content: 'e^x',
        expression: (val1, val2) => Math.exp(val1)
    },
    {
        value: 'custom',
        content: '2*e^x',
        expression: (val1, val2) => 2 * Math.exp(val1)
    },
    {
        value: 'tanh',
        content: 'tanh',
        expression: (val1, val2) => (Math.exp(2 * val1) - 1) / (Math.exp(2 * val1) + 1)
    },
]
