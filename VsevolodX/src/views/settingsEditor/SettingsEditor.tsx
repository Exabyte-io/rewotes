import React from 'react';
import { useSettings } from '../../context/SettingsContext';
import { Button, Card, Icon, Menu, MenuItem } from '@blueprintjs/core';
import styles from './SettingsEditor.module.scss';

// TODO: this is temporary file
const atomColors = [{'element': 'Si', 'color': 'blue'}, {'element': 'O', 'color': 'red'}, {'element': 'C', 'color': 'grey'}, {'element': 'H', 'color': 'black'}, {'element': 'Fe', 'color': 'orange'}]

function SettingsEditor() {
  const { settings, updateSettings } = useSettings();

  function toggleTheme() {
    const newTheme = settings.theme === 'light' ? 'dark' : 'light';
    updateSettings({ theme: newTheme });
  }

  return (
    <Card className={styles.SettingsEditor}>
      <h3>Settings</h3>
      <Button onClick={toggleTheme}>
        {settings.theme === 'light'? 'Dark' : 'Light'} Theme
      </Button>
      {/*TODO: Add ability to change colors in Settings */}
      <Menu>
        {atomColors.map((element) => {
          return(
            <MenuItem 
            text={element.element} 
            icon={<Icon icon='full-circle' color={element.color} />}
            />
          )
        })}
      </Menu>
    </Card>
  );
}

export default SettingsEditor;
