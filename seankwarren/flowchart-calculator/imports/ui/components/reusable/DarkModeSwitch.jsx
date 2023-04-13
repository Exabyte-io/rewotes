import React from 'react';
import Switch from 'react-switch';

const DarkModeSwitch = ({ isDarkMode, toggleDarkMode }) => {
    return (
        <Switch
            checked={isDarkMode}
            onChange={toggleDarkMode}
            offColor='#f0f0f0'
            onColor='#606060'
            offHandleColor='#ccc'
            onHandleColor='#444'
            handleDiameter={30}
            uncheckedIcon={false}
            checkedIcon={false}
            boxShadow='0px 1px 5px rgba(0, 0, 0, 0.6)'
            height={15}
            width={40}
            className='react-switch darkmode-switch'
        />
    );
};

export default DarkModeSwitch;
