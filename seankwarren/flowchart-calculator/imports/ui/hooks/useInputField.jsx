import React, { useState } from 'react'

const useInputField = () => {
    const [flowName, setFlowName] = useState('');

    const updateFlowName = (e) => {
        setFlowName(e.target.value);
    };

    return {flowName, updateFlowName}
}

export default useInputField