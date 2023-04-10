
import { customAlphabet } from 'nanoid'

const nanoid = customAlphabet('0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ', 6)

const createNode = (type, position, value, onInputChange) => {
    const id = `${type}${nanoid()}`;

    const newNode = {
        id,
        type: `${type}Node`,
        data: {
            value: value !== undefined ? value : getDefaultNodeValue(type),
            onChange: (newValue) => { 
                onInputChange(id, newValue) 
            },
        },
        position: position,
    };

    return newNode;
};

const getDefaultNodeValue = (type) => {
    switch (type) {
        case 'input':
            return 0;
        case 'binary':
            return 'add';
        case 'unary':
            return 'sin';
        case 'comparison':
            return 'greater';
        default:
            return null;
    }
};

export default createNode;