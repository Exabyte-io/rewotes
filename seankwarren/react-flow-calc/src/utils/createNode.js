import { nanoid } from 'nanoid';

const createNode = (type, position, value) => {
    const id = `${type}-${nanoid(11)}`;

    const newNode = {
        id,
        type: `${type}Node`,
        data: {
            value: value !== undefined ? value : getDefaultNodeValue(type),
            onChange: (value) => {
                // This function will be updated later in the component
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
        case 'comparison':
            return 'greater';
        default:
            return null;
    }
};

export default createNode;