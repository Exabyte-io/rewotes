import React from 'react';
import { useSettings } from '../../context/SettingsContext';
import { Button, Card, Icon, MenuItem } from '@blueprintjs/core';
import styles from './SettingsEditor.module.scss';
import Menu from '../../components/menu/Menu';
import VStack from '../../components/utils/VStack';
import ViewHeading from '../../components/view_heading/ViewHeading';
// TODO: this is temporary file

function SettingsEditor() {
  const { settings, updateSettings } = useSettings();
  const atomsData = useSettings().settings.atomsDisplayData || [];

  function toggleTheme() {
    const newTheme = settings.theme === 'light' ? 'dark' : 'light';
    updateSettings({ theme: newTheme });
  }

  async function loadAtomsData() {
    const res = await fetch('/atomsDisplayData.json');
    const data = await res.json();
    console.log(data);
    updateSettings({atomsDisplayData: data.atomsDisplayData});
  }

  return (
    <Card className={styles.SettingsEditor}>
      <VStack>
      <ViewHeading>
      <h4>Settings</h4>
      </ViewHeading>

      <Button 
        onClick={toggleTheme}
        icon={settings.theme === 'light'? 'moon' : 'lightbulb'}
      >
        {settings.theme === 'light'? 'Dark' : 'Light'} Theme
      </Button>
      {/*TODO: Add ability to change colors in Settings */}
      <Button onClick={() => loadAtomsData()} >Load atoms</Button>
      <Menu>
        {atomsData.map((element) => {
          return(
            <MenuItem 
            key={element.element}
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
