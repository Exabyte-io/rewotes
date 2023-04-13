import React, { useState } from 'react'

const useResizeable = () => {
    const [sizes, setSizes] = useState(['60%', '50%']);

    const handleSizeChange = (sizes) => {
        setSizes(sizes);
    };

    return {sizes, handleSizeChange}
}

export default useResizeable