import React from 'react';
import { useSettings } from '../../context/SettingsContext';
import { Button, Card, ControlGroup, Icon, MenuItem } from '@blueprintjs/core';
import styles from './SettingsEditor.module.scss';
import Menu from '../../components/menu/Menu';
import VStack from '../../components/utils/VStack';
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
      <h4>Settings</h4>
      <VStack>

      <Button 
        onClick={toggleTheme}
        icon={settings.theme === 'light'? 'moon' : 'lightbulb'}
      >
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
      </VStack>
    </Card>
  );
}

export default SettingsEditor;
