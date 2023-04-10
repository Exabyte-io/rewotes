import { nanoid } from 'nanoid';

const generateID = () => {
  const alphabet = '0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ';
  return nanoid(10, alphabet);
};

const createNode = (type, position, value, onInputChange) => {
    const id = `${type}${generateID()}`;

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