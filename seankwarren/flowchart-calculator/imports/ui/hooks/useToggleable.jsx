import React, { useState } from 'react'

const useToggleable = (val) => {
    const [isDarkMode, setIsDarkMode] = useState(!(!val))

    function toggleDarkMode() {
        setIsDarkMode(!isDarkMode);
    };

    return {isDarkMode, toggleDarkMode}
}

export default useToggleable